#pragma once
#include"include/SitranoHeader.h"
#include<vector>
#include<complex>
#include<fftw3.h>

class WTransform
{
public:
    WTransform(std::vector<float> input, float sampleRate)
    : signal(input), fs(sampleRate)
    {
        nfft = Sitrano::findNextPowerOfTwo(signal.size());
        signal.resize(nfft, 0.f);
        initFFTW();
    }

    virtual ~WTransform()
    {
        fftwf_destroy_plan(fwdPlan);
        fftwf_destroy_plan(invPlan);
        fftwf_free(inBuffer);
        fftwf_free(fDomainBuffer);
        fftwf_free(invBuffer);
    };

    virtual std::vector<float> process(const std::vector<float>& frequencies) = 0;

protected:
    float* inBuffer;
    fftwf_complex* fDomainBuffer;
    fftwf_complex* invBuffer;
    fftwf_plan fwdPlan;
    fftwf_plan invPlan;

    std::vector<float> signal;
    int nfft = 0;
    const float fs;

    void initFFTW()
    {
        inBuffer = (float*)fftwf_malloc(sizeof (float) * nfft);
        fDomainBuffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));
        invBuffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nfft);
        fwdPlan = fftwf_plan_dft_r2c_1d(nfft, inBuffer, fDomainBuffer, FFTW_MEASURE);
        invPlan = fftwf_plan_dft_1d(nfft, invBuffer, invBuffer, FFTW_BACKWARD, FFTW_MEASURE);
    }
    
    virtual std::vector<float> wavelet(float freq) = 0;
};