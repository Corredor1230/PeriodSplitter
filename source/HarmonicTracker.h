#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>
#include<iostream>
#include"SitranoHeader.h"

constexpr double M_PI = 3.14159265358979323846;

class HarmonicTracker {
public:

    HarmonicTracker(const Sitrano::AnalysisUnit& a,
        Sitrano::Results& r,
        std::vector<Sitrano::Peak> top,
        bool applyHann = true,
        float toleranceVal = 300.0);

    ~HarmonicTracker() {
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    }

    void analyze();
    //void analyzeFullFile();
    double interpolatePeak(int k, const std::vector<double>& mags);

    const std::vector<std::vector<float>>& getAmplitudes() const { return results.amps; }
    const std::vector<std::vector<float>>& getPhases() const { return results.phases; }
    const std::vector<std::vector<float>>& getFrequencies() const { return results.freqs; }

private:

    std::vector<float> window;
    float tolerance;        // Tolerance value used to find nearby bins
    bool applyHannWindow;

    fftwf_plan plan;
    float* input;
    std::vector<float> checker;
    fftwf_complex* output;

    std::vector<Sitrano::Peak> tFreqs;
    const Sitrano::AnalysisUnit& unit;
    Sitrano::Results& results;

    void initFFTW();
    void applyHann(float* data, int size);
};
