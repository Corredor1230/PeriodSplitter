#pragma once
#include"include/SitranoHeader.h"
#include<vector>
#include<fftw3.h>

class WTransform
{
public:
    WTransform(int windowSize);
    ~WTransform();

private:
    const int wSize;
    float* input;
    fftwf_complex* output;
    fftwf_plan plan;

    std::vector<float> wave(std::vector<float>& input, int windowSize);

    void initFFTW();

};