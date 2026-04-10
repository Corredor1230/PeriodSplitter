#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
#include"fftw3.h"
#include"include/SitranoHeader.h"

class MedianFilter{
public:
    MedianFilter(const int n_fft, const int hop);
    ~MedianFilter();

    struct HPSpectrogram{
        HPSpectrogram(Sitrano::Spectrogram h, Sitrano::Spectrogram p) : harmonic(h), percussive(p) {};
        Sitrano::Spectrogram harmonic;
        Sitrano::Spectrogram percussive;
    };

    struct ComplexHPSpec{
        ComplexHPSpec(Sitrano::ComplexSpectrogram h, Sitrano::ComplexSpectrogram p) : harmonic(h), percussive(p){};
        Sitrano::ComplexSpectrogram harmonic;
        Sitrano::ComplexSpectrogram percussive;
    };

    HPSpectrogram filter(const Sitrano::ComplexSpectrogram& input, const int filtSize);
    std::vector<float> processAudio(const std::vector<float>& input, const int filtSize);


private:

    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan fftPlan;

    float* ioutBuffer;
    fftwf_complex* invBuffer;
    fftwf_plan ifftPlan;

    const int nfft;
    const int hopSize;

    Sitrano::ComplexSpectrogram getComplexSpectrogram(const std::vector<float>& input, const int nfft, const int hopSize);

    ComplexHPSpec applyMask(HPSpectrogram& hp, Sitrano::ComplexSpectrogram& complex);

    std::vector<float> reconstructAudio(const Sitrano::ComplexSpectrogram& spec, const int nfft, const int hopSize);

    void initfftw(const int nfft);

};