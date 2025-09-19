#include"HarmonicTracker.h"

HarmonicTracker::HarmonicTracker(const std::vector<float>& file,
    const std::vector<int>& periodStarts, float sampleRate, float inPitch, 
    int numHarmonics, int fftSize, bool applyHann, float toleranceVal) :
    buffer(file),
    periods(periodStarts),
    pitch(inPitch),
    fs(sampleRate),
    N(fftSize),
    numHarmonics(numHarmonics),
    applyHannWindow(applyHann),
    tolerance(toleranceVal)
{
    window.resize(N);
    checker.resize(N * 2);
    for (int i = 0; i < N; i++)
    {
        window[i] = 0.0;
    }
    initFFTW();
}

void HarmonicTracker::initFFTW() {
    input = (float*)fftwf_malloc(sizeof(float) * N);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (N / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(N, input, output, FFTW_MEASURE);
}

void HarmonicTracker::applyHann(float* data, int size) {
    for (int i = 0; i < size; ++i) {
        float w = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (size - 1)));
        data[i] *= w;
    }
}

void HarmonicTracker::analyze() {
    size_t numFrames = periods.size() - 1;
    float ampThresh = dbToAmp(-50.0);
    harmonicAmplitudes.assign(numHarmonics, std::vector<float>(numFrames, 0.f));
    harmonicPhases.assign(numHarmonics, std::vector<float>(numFrames, 0.f));
    harmonicFrequencies.assign(numHarmonics, std::vector<float>(numFrames, 0.f));
    int startForTop = (int)((double)buffer.size() * 0.2);

    for (int j = 0; j < N * 2; ++j)
    {
        checker[j] = ((j + startForTop) < buffer.size()) ? buffer[j + startForTop] : 0.f;
    }

    std::vector<Peak> topFrequencies = findTopSorted(checker, N, fs, startForTop, true, tolerance);

    for (size_t i = 0; i < numFrames - 2; ++i) {
        int start = periods[i];
        int end = periods[i] + N;
        int periodLength = end - start;

        if (periodLength <= 0 || start + periodLength > buffer.size()) continue;

        // Repeat the period to fill the FFT input buffer
        for (int j = 0; j < N; ++j) {
            int srcIdx = start + (j % periodLength);
            input[j] = (srcIdx < buffer.size()) ? buffer[srcIdx] : 0.f;
        }

        if (applyHannWindow)
            applyHann(input, N);

        fftwf_execute(plan);

        // Approximate fundamental frequency for this period
        float f0 = float(fs) / float(periodLength);

        for (int h = 1; h <= topFrequencies.size(); ++h) {
            if (topFrequencies[h - 1].freq <= 20.0) continue;
            float targetFreq = topFrequencies[h - 1].freq;
            float binFreq = 0.f;
            if (i == 0) binFreq = targetFreq;
            else binFreq = findPeakWithinTolerance(harmonicFrequencies[h - 1][i - 1], tolerance, N, fs, output);
            
            int bin = static_cast<int>(std::round(targetFreq * N / float(fs)));

            if (bin < N / 2 + 1) {
                float real = output[bin][0];
                float imag = output[bin][1];

                float mag = std::sqrt(real * real + imag * imag);
                float amp = mag_to_amp(mag, N);
                float phase = std::atan2(imag, real);  // radians, range [-?, ?]

                harmonicAmplitudes[h - 1][i] = amp;
                harmonicPhases[h - 1][i] = phase;
                harmonicFrequencies[h - 1][i] = binFreq;
            }
        }
    }
}

double HarmonicTracker::interpolatePeak(int k, const std::vector<double>& mags)
{
    if (k <= 0 || k >= (int)mags.size() - 1) return 0.0;
    double m1 = mags[k - 1], m0 = mags[k], p1 = mags[k + 1];
    double denom = (m1 - 2 * m0 + p1);
    if (fabs(denom) < 1e-12) return 0.0;
    return 0.5 * (m1 - p1) / denom;
}

std::vector<HarmonicTracker::Peak> HarmonicTracker::findTopFrequencies(const std::vector<float>& input,
    int N, double sr)
{
    int N_out = N / 2 + 1;
    std::vector<fftw_complex> out(N_out);
    std::vector<double> in(N);

    // Hann window
    for (int i = 0; i < N; i++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
        in[i] = input[i] * w;
    }

    // FFTW plan
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in.data(), out.data(), FFTW_ESTIMATE);
    fftw_execute(plan);

    // Magnitudes
    std::vector<double> mags(N_out);
    for (int k = 0; k < N_out; k++) {
        double re = out[k][0];
        double im = out[k][1];
        mags[k] = sqrt(re * re + im * im);
    }

    // Collect peaks
    std::vector<Peak> peaks;
    peaks.reserve(N_out);
    for (int k = 1; k < N_out - 1; k++) { // skip DC and Nyquist
        double delta = interpolatePeak(k, mags);
        double freq = (k + delta) * sr / N;
        peaks.push_back({ freq, mags[k] });
    }

    // Sort by magnitude, take top 32
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
        return a.mag > b.mag;
        });
    if (peaks.size() > 32) peaks.resize(32);

    // Resort by frequency
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
        return a.freq < b.freq;
        });

    fftw_destroy_plan(plan);
    return peaks;
}

std::vector<HarmonicTracker::Peak> HarmonicTracker::findTopSorted(const std::vector<float>& input, int N, double fs,
    int startPos, bool useCentsTolerance, double tolValue, bool sumAmplitudesInCluster)
{
    const int outHarm = numHarmonics;
    const int Nout = N / 2 + 1;
    std::vector<double> in(N);
    int startSample = startPos;
    for (int n = 0; n < N; ++n) {
        double w = 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (N - 1)));
        in[n] = (n) < input.size() ? input[n] * w : 0.0;
    }

    std::vector<fftw_complex> out(Nout);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in.data(), out.data(), FFTW_ESTIMATE);
    fftw_execute(plan);

    std::vector<double> mags(Nout);
    for (int k = 0; k < Nout; ++k) {
        double re = out[k][0], im = out[k][1];
        mags[k] = std::sqrt(re * re + im * im);
    }

    std::vector<Peak> peaks;
    peaks.reserve(Nout);
    for (int k = 1; k < Nout - 1; ++k) {
        double delta = interp_delta(k, mags);
        double f = (k + delta) * fs / double(N);
        double amp = mag_to_amp(mags[k], N);
        peaks.push_back({ f, mags[k], amp });
    }

    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
        return a.amp > b.amp;
        });

    /*for (int i = 0; i < peaks.size(); i++)
    {
        if (i == 0 || i == peaks.size() - 1) continue;
        if (peaks[i].amp > 9.0e-5 && peaks[i].amp > peaks[i + 1].amp
            && peaks[i].amp > peaks[i - 1].amp) std::cout << peaks[i].amp << " F: " << peaks[i].freq << '\n';
    }*/

    std::vector<Cluster> clusters;
    clusters.reserve(Nout);
    auto withinTolerance = [&](double base, double comp) {
        double tolHz = cents_to_hz(base, tolValue);
        double diff = std::abs(comp - base);
        bool within = diff <= tolHz;
        return within;
        };

    /*for (int pFreq = 0; pFreq < peaks.size(); pFreq++)
    {
        float baseF = peaks[pFreq].freq;
        Cluster c;
        clusters.push_back(c);
        for (int cFreq = pFreq; cFreq < peaks.size(); cFreq++)
        {
            float compF = peaks[cFreq].freq;
            if (withinTolerance(baseF, compF))
            {
                clusters.back().add(peaks[cFreq]);
            }
            else
            {
                pFreq = cFreq;
                break;
            }
        }
    }

    std::sort(clusters.begin(), clusters.end(), [](const Cluster& a, const Cluster& b) {
        return a.maxAmp > b.maxAmp;
        });

    if (clusters.size() > outHarm) clusters.resize(outHarm);*/

    std::vector<Peak> merged;
    merged.reserve(outHarm);

    for (int i = 0; i < peaks.size() && merged.size() < outHarm; i++) {
        if (i == 0) {
            merged.push_back(peaks[i]);
        }
        else {
            bool newHarmonic = false;
            std::vector<Peak> temp = merged;
            for (auto& m : temp) {
                if (withinTolerance(peaks[i].freq, m.freq)) {
                    newHarmonic = false;
                    break;
                }
                else {
                    newHarmonic = true;
                    continue;
                }
            }
            if (newHarmonic && (peaks[i].freq > pitch * 0.8)
                && (peaks[i].freq < 20000.f)
                && (peaks[i].amp > 2.0e-6)) merged.push_back(peaks[i]);
        }
    }

    /*for (const auto& c : clusters) {
        Peak p = c.to_peak();
        if (!sumAmplitudesInCluster) p.amp = c.maxAmp;
        merged.push_back(p);
    }*/

    std::sort(merged.begin(), merged.end(), [](const Peak& a, const Peak& b) {
        return a.amp > b.amp;
        });

    if (merged.size() > outHarm) merged.resize(outHarm);

    // finally sort by frequency
    std::sort(merged.begin(), merged.end(), [](const Peak& a, const Peak& b) {
        return a.freq < b.freq;
        });

    fftw_destroy_plan(plan);
    return merged;
}

//void HarmonicTracker::analyzeFullFile()
//{
//    for (int period = 0; period < periods.size() - 1; period++)
//    {
//        int periodSize = periods[period + 1] - periods[period];
//        int pInW = std::floor<int>((float)N / periodSize);
//        std::vector<float> single(periodSize);
//
//        for (int i = 0; i < pInW; i++)
//        {
//            for (int samp = 0; samp < N; samp++)
//            {
//
//            }
//        }
//    }
//}

void HarmonicTracker::setNumHarmonics(int maxNum, bool adjustToPitch, float pitch)
{
    if (!adjustToPitch) numHarmonics = maxNum;
    else
    {
        int newMax = 0;
        for (int i = 0; i < maxNum; i++)
        {
            float highestHz = pitch * (i + 1);
            if (highestHz > fs / 2.0) break;
            newMax = i + 1;
        }
        numHarmonics = newMax;
    }
}


void HarmonicTracker::filterVector(std::vector<float>& input, int filterSize)
{
    std::vector<float> filter;
    filter.resize(filterSize);

    for (int i = 0; i < filter.size(); i++)
    {
        filter[i] = 0.f;
    }

    for (int i = 0; i < input.size(); i++)
    {
        int index = i % filterSize;
        filter[index] = input[i];
        float sum = 0.f;
        for (int j = 0; j < filter.size(); j++)
        {
            sum += filter[j];
        }
        if (filter.size() != 0) sum /= (float)filter.size();
        input[i] = sum;
    }
}