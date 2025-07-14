#include"HarmonicTracker.h"

HarmonicTracker::HarmonicTracker(const std::vector<float>& audioBuffer,
    const std::vector<int>& periodStarts, float sampleRate, int numHarmonics,
    int fftSize, bool applyHann) :
    buffer(audioBuffer),
    periods(periodStarts),
    fs(sampleRate),
    N(fftSize),
    numHarmonics(numHarmonics),
    applyHannWindow(applyHann)
{
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
    harmonicAmplitudes.assign(numHarmonics, std::vector<float>(numFrames, 0.f));
    harmonicPhases.assign(numHarmonics, std::vector<float>(numFrames, 0.f));

    for (size_t i = 0; i < numFrames; ++i) {
        int start = periods[i];
        int end = periods[i + 1];
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

        for (int h = 1; h <= numHarmonics; ++h) {
            float targetFreq = h * f0;
            int bin = static_cast<int>(std::round(targetFreq * N / float(fs)));

            if (bin < N / 2 + 1) {
                float real = output[bin][0];
                float imag = output[bin][1];

                float mag = std::sqrt(real * real + imag * imag);
                float phase = std::atan2(imag, real);  // radians, range [-?, ?]

                harmonicAmplitudes[h - 1][i] = mag;
                harmonicPhases[h - 1][i] = phase;
            }
        }
    }
}

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
    }
}