#include"SitranoAnalysis.h"

Analyzer::Analyzer(const Sitrano::AnalysisUnit& s, 
    Sitrano::Results& r, bool applyHann, float toleranceVal) :
    applyHannWindow(applyHann),
    tolerance(toleranceVal),
    ana(s),
    results(r)
{
    window.resize(ana.nfft);
    checker.resize(ana.nfft * 2);
    for (int i = 0; i < ana.nfft; i++)
    {
        window[i] = 0.0;
    }
    initFFTW();
}

void Analyzer::initFFTW() {
    input = (float*)fftwf_malloc(sizeof(float) * ana.nfft);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (ana.nfft / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(ana.nfft, input, output, FFTW_MEASURE);
}

void Analyzer::applyHann(double* data, int size) {
    for (int i = 0; i < size; ++i) {
        float w = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (size - 1)));
        data[i] *= w;
    }
}

void Analyzer::analyze(const Sitrano::Settings& settings) {

    int provisionalHop = ana.sampleRate / 60;

    if (settings.pitchAnalysis)
    {
        PitchFinder finder{ ana };
        results.pitch = finder.findPitch();
    }

    if (settings.periodAnalysis)
    {
        PeriodCutter cutter{ ana, results };
        if (results.pitch > 0.0)
            cutter.initialize(results.pitch);
        else
            cutter.initialize(ana.sampleRate / provisionalHop);

        results.sampleList = cutter.getCorrelationZeroes();
    }

    size_t numFrames = results.sampleList.size() > 1 ? results.sampleList.size() - 1 : provisionalHop;

    results.amps.assign(ana.numHarmonics, std::vector<float>(numFrames, 0.f));
    results.phases.assign(ana.numHarmonics, std::vector<float>(numFrames, 0.f));
    results.freqs.assign(ana.numHarmonics, std::vector<float>(numFrames, 0.f));

    int startForTop = (int)((double)ana.soundFile.size() * 0.12);

    for (int j = 0; j < ana.nfft * 2; ++j)
    {
        checker[j] = ((j + startForTop) < ana.soundFile.size()) ? ana.soundFile[j + startForTop] : 0.f;
    }

    std::vector<Sitrano::Peak> topFrequencies = findTopSorted(checker, ana.nfft, ana.sampleRate, startForTop, true, tolerance);

    if (settings.harmonicAnalysis)
    {
        HarmonicTracker tracker(
            ana,
            results,
            topFrequencies);

        tracker.analyze();
    }

}

std::vector<Sitrano::Peak> Analyzer::findTopSorted(const std::vector<double>& input, int N, double fs,
    int startPos, bool useCentsTolerance, double tolValue, bool sumAmplitudesInCluster)
{
    const int outHarm = ana.numHarmonics;
    const int Nout = N / 2 + 1;
    std::vector<double> in = input;
    int startSample = startPos;
    applyHann(in.data(), N);

    std::vector<fftw_complex> out(Nout);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in.data(), out.data(), FFTW_ESTIMATE);
    fftw_execute(plan);

    std::vector<double> mags(Nout);
    for (int k = 0; k < Nout; ++k) {
        double re = out[k][0], im = out[k][1];
        mags[k] = std::sqrt(re * re + im * im);
    }

    std::vector<Sitrano::Peak> peaks;
    peaks.reserve(Nout);
    for (int k = 1; k < Nout - 1; ++k) {
        double delta = Sitrano::interp_delta(k, mags);
        double f = (k + delta) * fs / double(N);
        double amp = Sitrano::mag_to_amp(mags[k], N);
        peaks.push_back({ f, mags[k], amp });
    }

    std::sort(peaks.begin(), peaks.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
        return a.amp > b.amp;
        });

    auto withinTolerance = [&](double base, double comp) {
        double tolHz = Sitrano::cents_to_hz(base, tolValue);
        double diff = std::abs(comp - base);
        bool within = diff <= tolHz;
        return within;
        };

    std::vector<Sitrano::Peak> merged;
    merged.reserve(outHarm);

    for (int i = 0; i < peaks.size() && merged.size() < outHarm; i++) {
        if (i == 0) {
            merged.push_back(peaks[i]);
        }
        else {
            bool newHarmonic = false;
            std::vector<Sitrano::Peak> temp = merged;
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
            if (newHarmonic && (peaks[i].freq > results.pitch * 0.9)
                && (peaks[i].freq < 20000.f)
                && (peaks[i].amp > 1.0e-6)) merged.push_back(peaks[i]);
        }
    }

    std::sort(merged.begin(), merged.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
        return a.amp > b.amp;
        });

    if (merged.size() > outHarm) merged.resize(outHarm);

    // finally sort by frequency
    std::sort(merged.begin(), merged.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
        return a.freq < b.freq;
        });

    fftw_destroy_plan(plan);
    return merged;
}

