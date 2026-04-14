#pragma once

#include "include/SitranoHeader.h"
#include <vector>
#include <fftw3.h>

class TransientAnalysis{
public:
    TransientAnalysis(const Sitrano::SingleTransientSettings& settings);
    ~TransientAnalysis();

private:
    void initfftw(int nfft);

    const Sitrano::SingleTransientSettings& settings;
    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan plan;
};