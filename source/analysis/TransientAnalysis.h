#pragma once

#include "include/SitranoHeader.h"
#include <vector>
#include <fftw3.h>

class TransientAnalysis{
public:
    TransientAnalysis(const Sihat::SingleTransientSettings& settings, const Sihat::AnalysisUnit& unit);

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
    std::vector<Sihat::OvertoneTrajectory> trackOvertonesInTime(const std::vector<float>& audio, const Sihat::SampleRange& range, const std::vector<Sihat::Peak>& targetPeaks, int numFrames);

    void initfftw(int nfft);

    const Sihat::SingleTransientSettings& settings;
    const Sihat::AnalysisUnit& u;
    float* inBuffer;
    fftwf_complex* outBuffer;
    fftwf_plan plan;
};