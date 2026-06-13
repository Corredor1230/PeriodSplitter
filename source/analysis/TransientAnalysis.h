#pragma once

#include "include/SitranoHeader.h"
#include "include/ResynthHeader.h"
#include <vector>
#include <fftw3.h>

class TransientAnalysis{
public:
    TransientAnalysis(const Sihat::SingleTransientSettings& settings, const Sihat::AnalysisUnit& unit, const float f0);

    Sihat::STransientResults analyze();

    ~TransientAnalysis();

private:

    Sihat::SampleRange findBoundaries(const std::vector<float>& input, const std::vector<Sihat::Sample>& env, const Sihat::SampleRange& range);
    Sihat::Spectrogram computeSpectrogram(const std::vector<float>& audio, const Sihat::SampleRange& range);
    std::vector<Sihat::Sample> getAmpEnvelope(const std::vector<float>& audio, const Sihat::SampleRange& range);
    int getRiseTime(const std::vector<Sihat::Sample>& env);
    std::vector<float> getSpectralCentroid(const Sihat::Spectrogram& spec);
    std::vector<float> getFlatnessCurve(const Sihat::Spectrogram& spec);
    std::vector<float> getBandEnvelopes(const Sihat::Spectrogram& spec);
    std::vector<Sihat::Peak> getMainOvertones(const std::vector<float>& audio, const Sihat::SampleRange& range);
    Sihat::TransientHarmonics trackOvertonesInTime(const std::vector<float>& audio, const Sihat::SampleRange& range, const std::vector<Sihat::Peak>& targetPeaks, int numFrames, bool trackFreq);
    std::vector<Synth::ModalComponent> extractModalComponents(const std::vector<float>& audio, int numFrames);

    void initfftw(int nfft);

    int findFirstPeak(const std::vector<Sihat::Sample>& env, int firstSample);

    const Sihat::SingleTransientSettings& settings;
    const Sihat::AnalysisUnit& u;
    const float f0;
    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan plan;
};