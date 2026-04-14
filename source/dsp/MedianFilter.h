#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
#include"fftw3.h"
#include"include/SitranoHeader.h"

class MedianFilter{
public:
    MedianFilter(Sitrano::HPSSSettings settings);
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

    HPSpectrogram filter(const Sitrano::ComplexSpectrogram& input);
    std::vector<float> processAudio(const std::vector<float>& input);


private:

    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan fftPlan;

    float* ioutBuffer;
    fftwf_complex* invBuffer;
    fftwf_plan ifftPlan;

    const int nfft;
    const int hopSize;
    int filtSize;

    Sitrano::ComplexSpectrogram getComplexSpectrogram(const std::vector<float>& input);

    ComplexHPSpec applyMask(HPSpectrogram& hp, Sitrano::ComplexSpectrogram& complex);

    std::vector<float> reconstructAudio(const Sitrano::ComplexSpectrogram& spec);

    void initfftw(const int nfft);

};