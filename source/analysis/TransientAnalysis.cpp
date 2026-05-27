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
    std::vector<Sihat::Sample> tEnv = getAmpEnvelope(u.soundFile, initRange);
    Sihat::SampleRange tRange = findBoundaries(u.soundFile, tEnv, initRange);
    float rms = Sihat::getRmsValue(u.soundFile, tRange.initSample, tRange.endSample);
    Sihat::Spectrogram tSpec = computeSpectrogram(u.soundFile, tRange);

    //Find representative frequencies
    std::vector<Sihat::Peak> stOvertones = getMainOvertones(u.soundFile, tRange);
    std::vector<float> stOvertonesAmps;
    for (int i = 0; i < stOvertones.size(); i++){stOvertonesAmps.push_back(stOvertones[i].amp);}

    //Track representative frequencies
    Sihat::TransientHarmonics tHarmonics = trackOvertonesInTime(u.soundFile, tRange, stOvertones, settings.overFrames);

    int tRise = getRiseTime(tEnv);
    std::vector<float> tCentroid = getSpectralCentroid(tSpec);
    std::vector<float> tFlatness = getFlatnessCurve(tSpec);
    std::vector<float> tBands = getBandEnvelopes(tSpec);

    int startFrame = tRange.initSample / settings.rmsHopSize;
    int endFrame = tRange.endSample / settings.rmsHopSize;
    startFrame = std::max(0, startFrame);
    endFrame = std::min(static_cast<int>(tEnv.size()), endFrame + 1);

    //Slicing envelopes to conform to actual transient
    std::vector<Sihat::Sample> slicedEnv(tEnv.begin() + startFrame, tEnv.begin() + endFrame);
    std::vector<float> slicedAEnv(slicedEnv.size());
    for (int i = 0; i < slicedEnv.size(); i++){slicedAEnv[i] = slicedEnv[i].value;}
    auto maxIt = std::max_element(slicedAEnv.begin(), slicedAEnv.end());
    float peakVal = *maxIt;

    // Sihat::STransientResults results {tRange, tRise, peakVal, slicedEnv, tCentroid, tFlatness, tBands, tHarmonics, rms, settings.rmsHopSize, settings.hopSize, tHarmonics.hopSize, settings.nfft, settings.nfft / 2 + 1};
    Sihat::STransientResults results{tRange, settings.rmsHopSize, settings.hopSize, tHarmonics.hopSize, settings.nfft, settings.nfft / 2 + 1, tRise, peakVal, slicedEnv, tCentroid, tFlatness, tBands, tHarmonics, rms};

    return results;
}

Sihat::SampleRange TransientAnalysis::findBoundaries(const std::vector<float>& input, const std::vector<Sihat::Sample>& tenv, const Sihat::SampleRange& range)
{

    std::vector<float> env(tenv.size());
    for (int i = 0; i < tenv.size(); i++){env[i] = tenv[i].value;}
    Sihat::SampleRange outRange;
    outRange.initSample = range.initSample;
    outRange.endSample = range.endSample;

    if (env.empty() || input.empty()) return outRange;

    auto maxIt = std::max_element(env.begin(), env.end());
    int peakFrame = static_cast<int>(std::distance(env.begin(), maxIt));
    float peakVal = *maxIt;

    float inThresh  = peakVal * std::pow(10.0f, settings.inThreshold / 20.0f);
    float outThresh = peakVal * std::pow(10.0f, settings.outThreshold / 20.0f);

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

    if (endFrame == -1) {
        endFrame = env.size() - 1; 
    }

    int appInitSamp = std::max(0, startFrame * settings.rmsHopSize);
    int appEndSamp = std::min(static_cast<int>(input.size() - 1), endFrame * settings.rmsHopSize);
    
    outRange.initSample = Sihat::findPreviousZero(input, appInitSamp);
    outRange.endSample = Sihat::findNextZero(input, appEndSamp);

    return outRange;
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

std::vector<Sihat::Sample> TransientAnalysis::getAmpEnvelope(const std::vector<float>& audio, const Sihat::SampleRange& range)
{
    int start = std::max(0, range.initSample);
    int end = std::min(static_cast<int>(audio.size()), range.endSample);
    int size = settings.rmsWindow;
    int hop = settings.rmsHopSize;
    std::vector<Sihat::Sample> env;

    for (int i = start; i < end; i += hop){
        float sum = 0.0f;
        int frameEnd = std::min(i + size, end);
        for (int j = i; j < frameEnd; j++)
        {
            sum += std::abs(audio[j]);
        }
        Sihat::Sample samp = {i, sum / (frameEnd - i)};
        env.push_back(samp);
    }

    return env;
}

int TransientAnalysis::getRiseTime(const std::vector<Sihat::Sample>& tenv)
{
    std::vector<float> env(tenv.size());
    for (int i = 0; i < tenv.size(); i++) {env[i] = tenv[i].value;}
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

std::vector<Sihat::Peak> TransientAnalysis::getMainOvertones(const std::vector<float>& audio, const Sihat::SampleRange& range)
{
    int maxOvertones = settings.maxOvertones;

    float* o_input;
    fftwf_complex* o_output;
    fftwf_plan o_plan;

    int nfft = Sihat::findNextPowerOfTwo(range.endSample - range.initSample);
    const int Nout = nfft / 2 + 1;

    o_input = (float*)fftwf_malloc(sizeof(float) * nfft);
    o_output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * Nout);
    
    // Note: If called frequently, consider making the plan a class member 
    // to avoid the expensive plan generation overhead.
    o_plan = fftwf_plan_dft_r2c_1d(nfft, o_input, o_output, FFTW_MEASURE);

    std::fill_n(o_input, nfft, 0.0f);
    std::copy(audio.begin() + range.initSample, audio.begin() + range.initSample + nfft, o_input);

    Sihat::applyHann<float>(o_input, nfft);

    fftwf_execute(o_plan);

    // 1. CALCULATE MAGNITUDES
    std::vector<double> mags(Nout);
    for (int k = 0; k < Nout; ++k) {
        double re = o_output[k][0];
        double im = o_output[k][1];
        mags[k] = std::sqrt(re * re + im * im);
    }

    // 2. CALCULATE GLOBAL NOISE FLOOR (Using all bins, not just peaks)
    // std::accumulate requires <numeric>
    double totalMagSum = std::accumulate(mags.begin(), mags.end(), 0.0);
    double averageFloor = totalMagSum / Nout;

    // 3. FIND & INTERPOLATE PEAKS 
    std::vector<Sihat::Peak> peaks;
    peaks.reserve(Nout); 
    for (int k = 1; k < Nout - 1; ++k) {
        // Find local maxima (basic peak-picking)
        if (mags[k] > mags[k - 1] && mags[k] > mags[k + 1])
        {
            double delta = Sihat::interp_delta(k, mags);
            double f = (k + delta) * u.sampleRate / double(nfft);
            double amp = Sihat::mag_to_amp(mags[k], nfft);
            double pha = std::atan2(o_output[k][1], o_output[k][0]);
            
            // Assuming your struct is { freq, mag, amp } based on previous usage
            peaks.push_back({ f, mags[k], amp, pha });
        }
    }

    // Sort all found peaks by magnitude (loudest first)
    std::sort(peaks.begin(), peaks.end(), [](const Sihat::Peak& a, const Sihat::Peak& b) {
        return a.mag > b.mag; 
    });

    std::vector<Sihat::Peak> merged;
    merged.reserve(settings.maxOvertones);

    std::vector<bool> peak_processed(peaks.size(), false);
    
    // Set a stricter crest factor for noisy signals
    const float minCrestFactor = 4.0f; 

    for (int i = 0; i < peaks.size() && merged.size() < settings.maxOvertones; ++i) {
        if (peak_processed[i]) {
            continue; // This peak was already added to a cluster
        }

        // This is the new "seed" peak
        const auto& seedPeak = peaks[i];
        float tolInHz = Sihat::cents_to_hz(seedPeak.freq, settings.tolInCents);
        peak_processed[i] = true;
        
        // Calculate crest factor against the GLOBAL average floor using Magnitudes
        float crestFactor = seedPeak.mag / averageFloor;

        // Apply filters to the seed peak
        if (seedPeak.freq < 500.0 ||
            seedPeak.freq > 20000.f ||
            crestFactor < minCrestFactor) {
            continue; // This peak is invalid (noise), skip it
        }

        bool withinTolerance = false;
        
        // Fixed: Changed loop index to 'j' to prevent shadowing the outer 'i'
        for (int j = 0; j < merged.size(); j++)
        {
            if (Sihat::withinTolerance(merged[j].freq, seedPeak.freq, tolInHz))
            {
                withinTolerance = true;
                break; // Small optimization: stop checking once we find a match
            }
        }

        if (!withinTolerance) {
            merged.push_back(seedPeak);
        }
    }

    // sort by frequency for a clean output
    std::sort(merged.begin(), merged.end(), [](const Sihat::Peak& a, const Sihat::Peak& b) {
        return a.freq < b.freq;
    });

    // Clean up FFTW resources
    fftwf_free(o_input);
    fftwf_free(o_output);
    fftwf_destroy_plan(o_plan);
    
    return merged;
}

Sihat::TransientHarmonics TransientAnalysis::trackOvertonesInTime(
    const std::vector<float>& audio, 
    const Sihat::SampleRange& range, 
    const std::vector<Sihat::Peak>& targetPeaks, 
    int numFrames)
{
    Sihat::TransientHarmonics tHarmonics;
    std::vector<float> floor;
    std::vector<Sihat::TransientOvertone> trajectories(targetPeaks.size());
    for (size_t i = 0; i < targetPeaks.size(); ++i) {
        trajectories[i].targetOvertone = targetPeaks[i];
        trajectories[i].envelope.reserve(numFrames); // Assuming 1 point per frame for cleaner structure
    }

    // Calculate our frame step (hop size)
    int totalRangeSamples = range.endSample - range.initSample;
    int frameStep = totalRangeSamples / numFrames;
    if (frameStep < 1) frameStep = 1;
    tHarmonics.hopSize = frameStep;
    int nfft = frameStep * 2;
    int sr = static_cast<int>(u.sampleRate);
    const int Nout = nfft / 2 + 1;

    float* input = (float*)fftwf_malloc(sizeof(float) * nfft);
    fftwf_complex* output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * Nout);
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(nfft, input, output, FFTW_MEASURE);

    std::vector<float> mags(Nout);


    // We start the window so the *center* of the first frame is at initSample
    int startOffset = range.initSample - (nfft / 2) > 0 ? range.initSample - (nfft / 2) : 0;
    tHarmonics.startSample = startOffset;

    const float trackingMinCrest = 1.5f; // Lower threshold to track the tail effectively

    for (int i = 0; i < numFrames; ++i) 
    {
        int winStart = startOffset + (i * frameStep);
        int winEnd = winStart + nfft;
        int frameSampleCenter = winStart + (nfft / 2); // Exact sample time for this frame

        // 1. Fill buffer with zero-padding logic
        std::fill_n(input, nfft, 0.0f);
        int copyStart = std::max(winStart, 0); // Assuming audio starts at 0
        int copyEnd = std::min(winEnd, (int)audio.size());

        if (copyStart < copyEnd) {
            int dstOffset = copyStart - winStart;
            std::copy(audio.begin() + copyStart, audio.begin() + copyEnd, input + dstOffset);
        }

        // Apply Hann Window using your exact logic
        for (int w = 0; w < nfft; ++w) {
            float winVal = 0.5f * (1.0f - std::cos(Sihat::TWO_PI * w / (nfft - 1)));
            input[w] *= winVal;
        }

        fftwf_execute(plan);

        // 2. Calculate Magnitudes and Noise Floor
        for (int k = 0; k < Nout; ++k) {
            double re = output[k][0];
            double im = output[k][1];
            mags[k] = std::sqrt(re * re + im * im);
        }

        //Calculate floor, but filter out main  harmonics
        float totalMagSum = std::accumulate(mags.begin(), mags.end(), 0.0);
        float averageFloor = totalMagSum / Nout;
        float ampFloor = Sihat::mag_to_amp(averageFloor, nfft);
        if (averageFloor < 1e-12) averageFloor = 1e-12; // Prevent div by zero
        floor.push_back(ampFloor);

        // 3. Analyze each target frequency
        for (size_t h = 0; h < targetPeaks.size(); ++h) 
        {
            float targetFreq = targetPeaks[h].freq;
            if (targetFreq <= 20.0f) continue; // Skip sub-audio

            // Determine search bounds based on tolerance
            float hzTolerance = Sihat::cents_to_hz(targetFreq, settings.o_tolInCents);
            int k_center = std::round(targetFreq * nfft / u.sampleRate);
            int k_spread = 2;
            
            // Constrain search bounds
            int k_min = std::max(1, k_center - k_spread);
            int k_max = std::min(Nout - 2, k_center + k_spread);

            int bestBin = -1;
            double maxMag = -1.0;

            // Search for the true peak within the tight tolerance window
            for (int k = k_min; k <= k_max; ++k) {
                if (mags[k] > mags[k - 1] && mags[k] > mags[k + 1]) {
                    if (mags[k] > maxMag) {
                        maxMag = mags[k];
                        bestBin = k;
                    }
                }
                else {
                    bestBin = Sihat::freqToBin(targetFreq, nfft, sr);
                    maxMag = mags[bestBin];
                }
            }

            Sihat::TrackedPoint point;

            // If we found a valid peak within tolerance
            if (bestBin != -1) {
                // Use your custom parabolic interpolation logic
                double m1 = mags[bestBin - 1], m0 = mags[bestBin], p1 = mags[bestBin + 1];
                double denom = (m1 - 2 * m0 + p1);
                double delta = 0.0;
                
                if (std::fabs(denom) >= 1e-12) {
                    delta = 0.5 * (m1 - p1) / denom;
                }

                point.freq = (bestBin + delta) * u.sampleRate / double(nfft);
                if (std::fabs(point.freq - targetFreq) > targetFreq * 0.01 && trajectories[h].envelope.size() > 1) {
                    point.freq = trajectories[h].envelope.back().freq;
                }
                point.crestFactor = maxMag / averageFloor;
                point.amp = Sihat::mag_to_amp(maxMag, Nout);
                
                // Determine if the peak is prominent enough to be considered "active"
                point.active = (point.crestFactor >= trackingMinCrest);
            } 
            else 
            {
                // Fallback: If no local max is found, write the target freq and 0 crest factor.
                // This prevents the frequency values from flying "all over the place" 
                // when the harmonic is missing.
                point.freq = targetFreq;
                point.crestFactor = 0.0;
                point.active = false;
            }

            trajectories[h].envelope.push_back(point);
        }
    }

    tHarmonics.overtones = trajectories;
    tHarmonics.floor = floor;

    fftwf_free(input);
    fftwf_free(output);
    fftwf_destroy_plan(plan);

    return tHarmonics;
}