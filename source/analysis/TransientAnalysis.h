#pragma once

#include "include/SitranoHeader.h"
#include <vector>
#include <fftw3.h>

class TransientAnalysis{
public:
    TransientAnalysis(const Sitrano::SingleTransientSettings& settings);
    ~TransientAnalysis();

private:

    Sitrano::SampleRange findBoundaries(const std::vector<float>& input);

    void initfftw(int nfft);

    const Sitrano::SingleTransientSettings& settings;
    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan plan;
};