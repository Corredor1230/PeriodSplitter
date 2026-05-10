#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
#include"fftw3.h"
#include"include/SitranoHeader.h"

class MedianFilter{
public:
    MedianFilter(Sihat::HPSSSettings settings);
    ~MedianFilter();

    struct HPSpectrogram{
        HPSpectrogram(Sihat::Spectrogram h, Sihat::Spectrogram p) : harmonic(h), percussive(p) {};
        Sihat::Spectrogram harmonic;
        Sihat::Spectrogram percussive;
    };

    struct ComplexHPSpec{
        ComplexHPSpec(Sihat::ComplexSpectrogram h, Sihat::ComplexSpectrogram p) : harmonic(h), percussive(p){};
        Sihat::ComplexSpectrogram harmonic;
        Sihat::ComplexSpectrogram percussive;
    };

    HPSpectrogram filter(const Sihat::ComplexSpectrogram& input);
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
    const float exponent;
    int filtSize;

    Sihat::ComplexSpectrogram getComplexSpectrogram(const std::vector<float>& input);

    ComplexHPSpec applyMask(HPSpectrogram& hp, Sihat::ComplexSpectrogram& complex);

    std::vector<float> reconstructAudio(const Sihat::ComplexSpectrogram& spec);

    void initfftw(const int nfft);

};