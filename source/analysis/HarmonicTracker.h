#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>
#include<iostream>
#include"include/SitranoHeader.h"
#include"support/IWindowStrategy.h"

class HarmonicTracker {
public:

    HarmonicTracker(const Sihat::AnalysisUnit& a,
        const Sihat::AnalysisConfig& conf,
        const std::vector<Sihat::Peak>& top,
        const std::vector<uint32_t>& sampleList,
        const int startSample);

    ~HarmonicTracker() {
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    }

    Sihat::HarmonicResults analyze();
    double interpolatePeak(int k, const std::vector<double>& mags);

    const std::vector<std::vector<float>>& getAmplitudes() const { return hResults.amps; }
    const std::vector<std::vector<float>>& getPhases() const { return hResults.phases; }
    const std::vector<std::vector<float>>& getFrequencies() const { return hResults.freqs; }

private:

    std::vector<float> window;
    std::unique_ptr<IWindowStrategy> mWindowStrategy;

    const int nfft;

    fftwf_plan plan;
    float* input;
    std::vector<float> checker;
    fftwf_complex* output;

    const std::vector<Sihat::Peak>& tFreqs;
    const Sihat::AnalysisUnit& unit;
    const Sihat::HarmonicSettings& settings;
    const Sihat::AnalysisConfig& config;
    const std::vector<uint32_t>& sList;
    const float sr;
    const int start;
    Sihat::HarmonicResults hResults;

    void initFFTW();
    void applyHann(float* data, int size);
    int findPeakSample();
};
