#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>
#include<iostream>
#include"HarmonicTracker.h"
#include"NoiseTracker.h"
#include"SitranoHeader.h"

constexpr double M_PI = 3.14159265358979323846;

class Analyzer {
public:

    Analyzer(Sitrano::AnalysisUnit,
        bool applyHann = true,
        float toleranceVal = 300.0);

    ~Analyzer() {
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    }

    void analyze();
    void setNumHarmonics(int maxNum, bool adjustToPitch, float pitch);

    const std::vector<std::vector<float>>& getAmplitudes() const { return harmonicAmplitudes; }
    const std::vector<std::vector<float>>& getPhases() const { return harmonicPhases; }
    const std::vector<std::vector<float>>& getFrequencies() const { return harmonicFrequencies; }

    void filterVector(std::vector<float>& input, int filterSize);

private:

    std::vector<Sitrano::Peak> findTopSorted(
        const std::vector<float>& input,
        int N, double fs, int startPos,
        bool useCentsTolerance = true,
        double tolValue = 15.0,
        bool sumAmplitudesInCluster = true);

    const std::vector<float>& buffer;
    const std::vector<int>& periods;
    std::vector<float> window;
    float fs;               // Sampling rate
    float tolerance;        // Tolerance value used to find nearby bins
    float pitch;
    int N;                // FFT size
    int numHarmonics;
    bool applyHannWindow;

    Sitrano::AnalysisUnit ana;
    Sitrano::Results results;

    fftwf_plan plan;
    float* input;
    std::vector<float> checker;
    fftwf_complex* output;

    std::vector<std::vector<float>> harmonicAmplitudes;
    std::vector<std::vector<float>> harmonicPhases;
    std::vector<std::vector<float>> harmonicFrequencies;

    void initFFTW();
    void applyHann(float* data, int size);
};
