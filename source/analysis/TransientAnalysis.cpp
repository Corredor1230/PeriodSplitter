#pragma once

#include"TransientAnalysis.h"
#include<algorithm>
#include<functional>

TransientAnalysis::TransientAnalysis(const Sitrano::SingleTransientSettings& settings) : settings(settings)
{
    initfftw(settings.nfft);
}


TransientAnalysis::~TransientAnalysis()
{
    fftwf_free(inBuffer);
    fftwf_free(outBuffer);
    fftwf_destroy_plan(plan);
}

void TransientAnalysis::initfftw(int nfft)
{
    inBuffer = (float*)fftwf_malloc(sizeof(float) * nfft);
    outBuffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));

    plan = fftwf_plan_dft_r2c_1d(nfft, inBuffer, outBuffer, FFTW_MEASURE);
}

Sitrano::SampleRange TransientAnalysis::findBoundaries(const std::vector<float>& input)
{
    Sitrano::SampleRange range;
    range.initSample = settings.startSample;

    std::vector<float> rmsList;
    bool transientFound = false;
    int tInitSample = 0;
    int tEndSample = 0;
    float rmsPeak = 0.f;
    int numBelowPeak = 0;

    std::vector<float> rmsWindow(settings.rmsWindow); // This will hold the transient

    for (int samp = settings.startSample; samp < (input.size() / 4); samp += settings.rmsHopSize)
    {
        float windowSum = 0.f;

        for (int rmsSamp = 0; rmsSamp < settings.rmsWindow && (rmsSamp + samp) < input.size(); rmsSamp++)
        {
            float x = input[samp + rmsSamp];
            windowSum += (x * x);
            rmsWindow[rmsSamp] = x;
        }

        float rms = std::sqrt(windowSum / settings.rmsWindow);
        if (rms > rmsPeak) rmsPeak = rms;
        rmsList.push_back(rms);
    }

    std::transform(rmsList.begin(), rmsList.end(), rmsList.begin(), [rmsPeak](float f){return f / rmsPeak;});

    for (int i = 0; i < rmsList.size(); i++)
    {
        float currVal = rmsList[i];
    }

    if (!transientFound) {
        // No transient, just find the first peak in a small window
        int peakSample = Sitrano::findPeakSample(input, 0, std::min((int)input.size(), 8192), true);
        range.endSample = Sitrano::findPreviousZero(input, peakSample);
        return range;
    }

    range.initSample = Sitrano::findPreviousZero(input, tInitSample);
    range.endSample = Sitrano::findPreviousZero(input, tEndSample);
    return range;
}