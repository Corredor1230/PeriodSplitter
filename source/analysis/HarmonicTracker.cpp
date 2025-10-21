#include"HarmonicTracker.h"

HarmonicTracker::HarmonicTracker(
    const Sitrano::AnalysisUnit& a,
    const Sitrano::HarmonicSettings& h,
    Sitrano::Results& r,
    std::vector<Sitrano::Peak> top) :
    unit(a),
    settings(h),
    results(r),
    tFreqs(top)
{
    window.resize(unit.nfft);
    checker.resize(unit.nfft * 2);
    for (int i = 0; i < unit.nfft; i++)
    {
        window[i] = 0.0;
    }
    initFFTW();
}

void HarmonicTracker::initFFTW() {
    input = (float*)fftwf_malloc(sizeof(float) * unit.nfft);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (unit.nfft / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(unit.nfft, input, output, FFTW_MEASURE);
}

void HarmonicTracker::applyHann(float* data, int size) {
    for (int i = 0; i < size; ++i) {
        float w = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (size - 1)));
        data[i] *= w;
    }
}

void HarmonicTracker::analyze() {
    size_t numFrames{0};
    int frameStep{ 1 };

    switch (settings.style)
    {
    case Sitrano::WindowStyle::periodLoop:
    {
        frameStep = 2;
        numFrames = (results.sampleList.size() - 1) / frameStep;
        break;
    }
    case Sitrano::WindowStyle::singlePeriod:
    {
        frameStep = 2;
        numFrames = (results.sampleList.size() - 1) / frameStep;
        break;
    }
    case Sitrano::WindowStyle::audioChunk:
    {
        frameStep = 1;
        int chunkSize = unit.soundFile.size() - results.sampleList[0];
        numFrames = chunkSize / unit.hopSize;
        break;
    }
    default:
    {
        frameStep = 2;
        numFrames = (results.sampleList.size() - 1) / frameStep;
        break;
    }
    }

    float ampThresh = Sitrano::dbToAmp(-50.0);

    for (size_t i = 0; i < numFrames - (frameStep + 1); i+=frameStep) {

        int start{ 0 };
        int end{ 0 };
        int periodLength{ 0 };

        //results.finalSamples.push_back(results.sampleList[i]);

        switch (settings.style)
        {
        case Sitrano::WindowStyle::periodLoop:
        {
            start = results.sampleList[i];
            end = results.sampleList[i + frameStep];
            periodLength = end - start;

            if (periodLength <= 0 || end > unit.soundFile.size()) continue;

            // Repeat the period to fill the FFT input buffer
            for (int j = 0; j < unit.nfft; ++j) {
                int srcIdx = start + (j % periodLength);
                input[j] = (srcIdx < unit.soundFile.size()) ? unit.soundFile[srcIdx] : 0.f;
            }
            break;
        }
        case Sitrano::WindowStyle::singlePeriod:
        {
            start = results.sampleList[i];
            end = results.sampleList[i + frameStep];
            periodLength = end - start;

            if (periodLength <= 0 || end > unit.soundFile.size()) continue;

            for (int j = 0; j < unit.nfft; ++j) {
                int srcIdx = start + j;
                input[j] = (srcIdx < end) ? unit.soundFile[srcIdx] : 0.f;
            }
            break;
        }
        case Sitrano::WindowStyle::audioChunk:
        {
            start = i * unit.hopSize;
            end = start + unit.nfft;
            results.finalSamples.push_back(start);

            if (end > unit.soundFile.size()) continue;

            for (int j = 0; j < unit.nfft; ++j) {
                input[j] = unit.soundFile[start + j];
            }
            break;
        }
        default:
        {
            start = results.sampleList[i];
            end = results.sampleList[i + frameStep];
            periodLength = end - start;

            if (periodLength <= 0 || end > unit.soundFile.size()) continue;

            // Repeat the period to fill the FFT input buffer
            for (int j = 0; j < unit.nfft; ++j) {
                int srcIdx = start + (j % periodLength);
                input[j] = (srcIdx < unit.soundFile.size()) ? unit.soundFile[srcIdx] : 0.f;
            }
            break;
        }
        }
        if (settings.applyHanning)
            applyHann(input, unit.nfft);

        fftwf_execute(plan);

        // Approximate fundamental frequency for this period
        float f0 = float(unit.sampleRate) / float(periodLength);
        if (settings.style != Sitrano::WindowStyle::audioChunk)
        {
            for (int b = 0; b < frameStep && i + b < results.sampleList.size(); b++)
            {
                results.finalSamples.push_back(results.sampleList[i + b]);
            }
        }
        
        for (int h = 1; h <= tFreqs.size(); ++h) {
            if (tFreqs[h - 1].freq <= 20.0) continue;
            float targetFreq = tFreqs[h - 1].freq;
            float binFreq = 0.f;
            binFreq = Sitrano::findPeakWithinTolerance(targetFreq, settings.toleranceValue, unit.nfft, unit.sampleRate, output);
            
            int bin = static_cast<int>(std::round(targetFreq * unit.nfft / float(unit.sampleRate)));

            if (bin < unit.nfft / 2 + 1) {
                float real = output[bin][0];
                float imag = output[bin][1];

                float mag = std::sqrt(real * real + imag * imag);
                float amp = Sitrano::mag_to_amp(mag, unit.nfft);
                float phase = std::atan2(imag, real);  // radians, range [-?, ?]
                
                for (int step = 0; step < frameStep; step++)
                {
                    results.amps[h - 1].push_back(amp);
                    results.phases[h - 1].push_back(phase);
                    results.freqs[h - 1].push_back(binFreq);
                    
                }
            }
        }
    }
    if (results.finalSamples.size() > results.amps[0].size())
    {
        results.finalSamples.resize(results.amps[0].size());
    }
    else if (results.finalSamples.size() < results.amps[0].size())
    {
        for (int i = results.finalSamples.size(); i < results.amps[0].size(); i++)
        {
            std::cerr << "This is wrong\n";
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