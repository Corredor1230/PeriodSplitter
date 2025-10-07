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
    }
    case Sitrano::WindowStyle::singlePeriod:
    {
        frameStep = 2;
        numFrames = (results.sampleList.size() - 1) / frameStep;
    }
    case Sitrano::WindowStyle::audioChunk:
    {
        frameStep = 1;
        int chunkSize = unit.soundFile.size() - results.sampleList[0];
        numFrames = chunkSize / unit.hopSize;
    }
    }

    float ampThresh = Sitrano::dbToAmp(-50.0);

    for (size_t i = 0; i < numFrames - 2; i+=frameStep) {

        int start{ 0 };
        int end{ 0 };
        int periodLength{ 0 };

        switch (settings.style)
        {
            case Sitrano::WindowStyle::periodLoop:
            {
                start = results.sampleList[i];
                end = results.sampleList[i + 3];
                periodLength = end - start;

                if (periodLength <= 0 || end > unit.soundFile.size()) continue;

                // Repeat the period to fill the FFT input buffer
                for (int j = 0; j < unit.nfft; ++j) {
                    int srcIdx = start + (j % periodLength);
                    input[j] = (srcIdx < unit.soundFile.size()) ? unit.soundFile[srcIdx] : 0.f;
                }
            }
            case Sitrano::WindowStyle::singlePeriod:
                start = results.sampleList[i];
                end = results.sampleList[i + 1];
        }
        if (settings.applyHanning)
            applyHann(input, unit.nfft);

        fftwf_execute(plan);

        // Approximate fundamental frequency for this period
        float f0 = float(unit.sampleRate) / float(periodLength);

        for (int h = 1; h <= tFreqs.size(); ++h) {
            if (tFreqs[h - 1].freq <= 20.0) continue;
            float targetFreq = tFreqs[h - 1].freq;
            float binFreq = 0.f;
            if (i == 0) binFreq = targetFreq;
            else binFreq = Sitrano::findPeakWithinTolerance(results.freqs[h - 1][i - 1], settings.toleranceValue, unit.nfft, unit.sampleRate, output);
            
            int bin = static_cast<int>(std::round(targetFreq * unit.nfft / float(unit.sampleRate)));

            if (bin < unit.nfft / 2 + 1) {
                float real = output[bin][0];
                float imag = output[bin][1];

                float mag = std::sqrt(real * real + imag * imag);
                float amp = Sitrano::mag_to_amp(mag, unit.nfft);
                float phase = std::atan2(imag, real);  // radians, range [-?, ?]

                results.amps[h - 1][i] = amp;
                results.phases[h - 1][i] = phase;
                results.freqs[h - 1][i] = binFreq;
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