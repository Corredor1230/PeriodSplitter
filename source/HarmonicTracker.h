#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>

constexpr double M_PI = 3.14159265358979323846;

class HarmonicTracker {
public:
    HarmonicTracker(const std::vector<float>& audioBuffer,
        const std::vector<int>& periodStarts,
        float sampleRate,
        int numHarmonics = 16,
        int fftSize = 4096,
        bool applyHann = true);

    ~HarmonicTracker() {
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    }

    void analyze();
    void setNumHarmonics(int maxNum, bool adjustToPitch, float pitch);

    const std::vector<std::vector<float>>& getAmplitudes() const { return harmonicAmplitudes; }
    const std::vector<std::vector<float>>& getPhases() const { return harmonicPhases; }

private:
    const std::vector<float>& buffer;
    const std::vector<int>& periods;
    float fs;               // Sampling rate
    int N;                // FFT size
    int numHarmonics;
    bool applyHannWindow;

    fftwf_plan plan;
    float* input;
    fftwf_complex* output;

    std::vector<std::vector<float>> harmonicAmplitudes;
    std::vector<std::vector<float>> harmonicPhases;

    void initFFTW();
    void applyHann(float* data, int size);
};
