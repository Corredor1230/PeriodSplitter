#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>
#include<iostream>
#include"include/SitranoHeader.h"
#include"support/IWindowStrategy.h"

//constexpr double M_PI = 3.14159265358979323846;

class HarmonicTracker {
public:

    HarmonicTracker(const Sitrano::AnalysisUnit& a,
        const Sitrano::AnalysisConfig& conf,
        const std::vector<Sitrano::Peak>& top,
        const std::vector<uint32_t>& sampleList);

    ~HarmonicTracker() {
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    }

    void analyze();
    Sitrano::HarmonicResults getEnvelopes();
    //void analyzeFullFile();
    double interpolatePeak(int k, const std::vector<double>& mags);

    const std::vector<std::vector<float>>& getAmplitudes() const { return hResults.amps; }
    const std::vector<std::vector<float>>& getPhases() const { return hResults.phases; }
    const std::vector<std::vector<float>>& getFrequencies() const { return hResults.freqs; }

private:

    std::vector<float> window;
    std::unique_ptr<IWindowStrategy> mWindowStrategy;

    int nfft;

    fftwf_plan plan;
    float* input;
    std::vector<float> checker;
    fftwf_complex* output;

    const std::vector<Sitrano::Peak>& tFreqs;
    const Sitrano::AnalysisUnit& unit;
    const Sitrano::HarmonicSettings& settings;
    const Sitrano::AnalysisConfig& config;
    const std::vector<uint32_t>& sList;
    Sitrano::HarmonicResults hResults;

    void initFFTW();
    void applyHann(float* data, int size);
};
