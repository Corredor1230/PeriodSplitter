#pragma once

#include"TransientAnalysis.h"
#include<algorithm>
#include<functional>

TransientAnalysis::TransientAnalysis(const Sihat::SingleTransientSettings& settings, const Sihat::AnalysisUnit& unit) : settings(settings), u(unit)
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

Sihat::STransientResults TransientAnalysis::analyze()
{
    Sihat::SampleRange initRange{0, u.soundFile.size() / 4};
    std::vector<float> tEnv = getAmpEnvelope(u.soundFile, initRange);
    Sihat::SampleRange tRange = findBoundaries(u.soundFile, tEnv, initRange);
    float rms = Sihat::getRmsValue(u.soundFile, tRange.initSample, tRange.endSample);
    Sihat::Spectrogram tSpec = computeSpectrogram(u.soundFile, tRange);
    int tRise = getRiseTime(tEnv);
    std::vector<float> tCentroid = getSpectralCentroid(tSpec);
    std::vector<float> tFlatness = getFlatnessCurve(tSpec);
    std::vector<float> tBands = getBandEnvelopes(tSpec);

    int startFrame = tRange.initSample / settings.rmsHopSize;
    int endFrame = tRange.endSample / settings.rmsHopSize;
    startFrame = std::max(0, startFrame);
    endFrame = std::min(static_cast<int>(tEnv.size()), endFrame + 1);

    std::vector<float> slicedEnv(tEnv.begin() + startFrame, tEnv.begin() + endFrame);

    auto maxIt = std::max_element(slicedEnv.begin(), slicedEnv.end());
    float peakVal = *maxIt;

    Sihat::STransientResults results {tRange, tRise, peakVal, slicedEnv, tCentroid, tFlatness, tBands, rms, settings.rmsHopSize, settings.hopSize, settings.nfft, settings.nfft / 2 + 1};

    return results;
}

Sihat::SampleRange TransientAnalysis::findBoundaries(const std::vector<float>& input, const std::vector<float>& env, const Sihat::SampleRange& range)
{

    Sihat::SampleRange outRange;
    outRange.initSample = range.initSample;
    outRange.endSample = range.endSample;

    if (env.empty() || input.empty()) return outRange;

    auto maxIt = std::max_element(env.begin(), env.end());
    int peakFrame = static_cast<int>(std::distance(env.begin(), maxIt));
    float peakVal = *maxIt;

    float inThresh = peakVal * settings.inThreshold;
    float outThresh = peakVal * settings.outThreshold;

    int startFrame = 0;
    for (int i = peakFrame; i >= 0; i--)
    {
        if (env[i] <= inThresh)
        {
            startFrame = i;
            break;
        }
    }

    int endFrame = -1;
    for (int i = peakFrame; i < env.size(); i++)
    {
        if (env[i] <= outThresh)
        {
            endFrame = i;
            break;
        }
    }

    int appInitSamp = std::max(0, startFrame * settings.rmsHopSize);
    int appEndSamp = std::min(static_cast<int>(input.size() - 1), endFrame * settings.rmsHopSize);
    
    outRange.initSample = Sihat::findPreviousZero(input, appInitSamp);
    outRange.endSample = Sihat::findNextZero(input, appEndSamp);

    return outRange;

    // std::vector<float> rmsList;
    // bool transientFound = false;
    // int tInitSample = 0;
    // int tEndSample = 0;
    // float rmsPeak = 0.f;
    // int numBelowPeak = 0;

    // std::vector<float> rmsWindow(settings.rmsWindow); // This will hold the transient
    // bool audioDetected = false;

    // for (int samp = settings.startSample; samp < (input.size() / 4); samp += settings.rmsHopSize)
    // {
    //     if (input[samp] == 0.f && !audioDetected) continue;
    //     else if (input[samp] != 0.f && !audioDetected) audioDetected = true;
        
    //     float windowSum = 0.f;

    //     for (int rmsSamp = 0; rmsSamp < settings.rmsWindow && (rmsSamp + samp) < input.size(); rmsSamp++)
    //     {
    //         float x = input[samp + rmsSamp];
    //         windowSum += (x * x);
    //         rmsWindow[rmsSamp] = x;
    //     }

    //     float rms = std::sqrt(windowSum / settings.rmsWindow);
    //     if (rms > rmsPeak) rmsPeak = rms;
    //     rmsList.push_back(rms);
    // }

    // std::transform(rmsList.begin(), rmsList.end(), rmsList.begin(), [rmsPeak](float f){return f / rmsPeak;});

    // for (int i = 0; i < rmsList.size(); i++)
    // {
    //     float currVal = rmsList[i];
    // }

    // if (!transientFound) {
    //     // No transient, just find the first peak in a small window
    //     int peakSample = Sitrano::findPeakSample(input, 0, std::min((int)input.size(), 8192), true);
    //     range.endSample = Sitrano::findPreviousZero(input, peakSample);
    //     return range;
    // }

    // range.initSample = Sitrano::findPreviousZero(input, tInitSample);
    // range.endSample = Sitrano::findPreviousZero(input, tEndSample);
    // return range;
}

Sihat::Spectrogram TransientAnalysis::computeSpectrogram(const std::vector<float>& audio, const Sihat::SampleRange& range){
    int start = std::max(0, range.initSample);
    int end = std::min(static_cast<int>(audio.size()), range.endSample);
    int length = end - start;
    int numBins = settings.nfft / 2 + 1;

    if(length < settings.nfft) {
        return {{}, 0, numBins};
    }

    std::vector<float> window = Sihat::getHannWindow(settings.nfft);
    int numFrames = 1 + (length - settings.nfft) / settings.hopSize;

    Sihat::Spectrogram spec;
    spec.numFrames = numFrames;
    spec.numBins = numBins;
    spec.data.resize(numFrames * numBins);

    for (int i = 0; i < numFrames; i++)
    {
        int frameStart = start + i * settings.hopSize;
        for (int j = 0; j < settings.nfft; j++)
        {
            inBuffer[j]  = audio[frameStart + j] * window[j];
        }

        fftwf_execute(plan);

        for (int j = 0; j < numBins; j++)
        {
            float re = outBuffer[j][0];
            float im = outBuffer[j][1];
            spec.at(i, j) = std::sqrt(re * re + im * im);
        }
    }

    return spec;
}

std::vector<float> TransientAnalysis::getAmpEnvelope(const std::vector<float>& audio, const Sihat::SampleRange& range)
{
    int start = std::max(0, range.initSample);
    int end = std::min(static_cast<int>(audio.size()), range.endSample);
    int size = settings.rmsWindow;
    int hop = settings.rmsHopSize;
    std::vector<float> env;

    for (int i = start; i < end; i += hop){
        float sum = 0.0f;
        int frameEnd = std::min(i + size, end);
        for (int j = i; j < frameEnd; j++)
        {
            sum += std::abs(audio[j]);
        }
        env.push_back(sum / (frameEnd - i));
    }

    return env;
}

int TransientAnalysis::getRiseTime(const std::vector<float>& env)
{
    if(env.empty()) return 0;

    auto maxIt = std::max_element(env.begin(), env.end());
    int maxIdx = static_cast<int>(std::distance(env.begin(), maxIt));
    int startIdx = 0;
    float threshold = (*maxIt) * 0.1f;
    for (int i = maxIdx; i >= 0; i--)
    {
        if (env[i] < threshold)
        {
            startIdx = i;
            break;
        }
    }
    return (maxIdx - startIdx) * settings.rmsHopSize;
}

std::vector<float> TransientAnalysis::getSpectralCentroid(const Sihat::Spectrogram& spec)
{
    std::vector<float> centroids(spec.numFrames);
    float nyquist = u.sampleRate / 2.0f;

    for (int i = 0; i < spec.numFrames; i++)
    {
        float num = 0.0f;
        float den = 0.0f;
        for (int j = 0; j < spec.numBins; j++)
        {
            float freq = (j * nyquist) / (spec.numBins - 1);
            float mag = spec.at(i, j);
            num += freq * mag;
            den += mag;
        }

        centroids[i] = (den > 1e-6f) ? (num / den) : 0.0f;
    }

    return centroids;
}

std::vector<float> TransientAnalysis::getFlatnessCurve(const Sihat::Spectrogram& spec)
{
    std::vector<float> flatness(spec.numFrames);

    for (int i = 0; i < spec.numFrames; i++)
    {
        float sumMag = 0.0f;
        float sumLogMag = 0.0f;
        int validBins = 0;

        for (int j = 0; j < spec.numBins; j++)
        {
            float mag = spec.at(i, j);
            if (mag > 1e-6f)
            {
                sumMag += mag;
                sumLogMag += std::log(mag);
                validBins++;
            }
        }

        if (validBins > 0 && sumMag > 1e-6f)
        {
            float geomMean = std::exp(sumLogMag / validBins);
            float arithMean = sumMag / validBins;
            flatness[i] = geomMean / arithMean;
        }
        else
        {
            flatness[i] = 0.0f;
        }
    }

    return flatness;
}

std::vector<float> TransientAnalysis::getBandEnvelopes(const Sihat::Spectrogram& spec)
{
    if (spec.numFrames == 0) return {};
    int numBands = settings.numBands;
    std::vector<float> envelopes(numBands * spec.numFrames, 0.0f);
    int binsPerBand = spec.numBins / numBands;

    for (int i = 0; i < spec.numFrames; i++)
    {
        for (int b = 0; b < numBands; b++)
        {
            float bandEnergy = 0.0f;
            int startBin = b * binsPerBand;
            int endBin = (b == numBands - 1) ? spec.numBins : (b + 1) * binsPerBand;

            for (int j = startBin; j < endBin; j++)
            {
                float mag = spec.at(i, j);
                bandEnergy += mag * mag;
            }
            envelopes[b * spec.numFrames + i] = bandEnergy;
        }
    }

    return envelopes;
}

