#pragma once

#include"TransientAnalysis.h"
#include<algorithm>
#include<functional>
#include<Eigen/Dense>

TransientAnalysis::TransientAnalysis(const Sihat::SingleTransientSettings& settings, const Sihat::AnalysisUnit& unit, const float pitch) : settings(settings), u(unit), f0(pitch)
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

    Synth::STransient str;

    Sihat::SampleRange initRange{0, u.soundFile.size() / 4};
    std::vector<Sihat::Sample> tEnv = getAmpEnvelope(u.soundFile, initRange);
    Sihat::SampleRange tRange = findBoundaries(u.soundFile, tEnv, initRange);
    int firstPeakSample = findFirstPeak(tEnv, tRange.initSample);
    int newinit = Sihat::findPreviousZero(u.soundFile, firstPeakSample);
    if (newinit < tRange.endSample) tRange.initSample = newinit;
    tEnv = getAmpEnvelope(u.soundFile, tRange);
    float rms = Sihat::getRmsValue(u.soundFile, tRange.initSample, tRange.endSample);
    Sihat::Spectrogram tSpec = computeSpectrogram(u.soundFile, tRange);

    //Find representative frequencies
    std::vector<Sihat::Peak> stOvertones = getMainOvertones(u.soundFile, tRange);
    std::vector<float> stOvertonesAmps;
    for (int i = 0; i < stOvertones.size(); i++){stOvertonesAmps.push_back(stOvertones[i].amp);}

    //Track representative frequencies
    Sihat::TransientHarmonics tHarmonics = trackOvertonesInTime(u.soundFile, tRange, stOvertones, settings.overFrames, false);

    int tRise = getRiseTime(tEnv);
    std::vector<float> tCentroid = getSpectralCentroid(tSpec);
    std::vector<float> tFlatness = getFlatnessCurve(tSpec);
    std::vector<float> tBands = getBandEnvelopes(tSpec);

    int startFrame = tRange.initSample / settings.rmsHopSize;
    int endFrame = tRange.endSample / settings.rmsHopSize;
    startFrame = std::max(0, startFrame);
    endFrame = std::min(static_cast<int>(tEnv.size()), endFrame + 1);

    //Slicing envelopes to conform to actual transient
    std::vector<float> tAEnv(tEnv.size());
    for (int i = 0; i < tEnv.size(); i++){tAEnv[i] = tEnv[i].value;}
    auto maxIt = std::max_element(tAEnv.begin(), tAEnv.end());
    float peakVal = *maxIt;

    tRise = firstPeakSample - tRange.initSample;
    std::vector<float> audioSlice = Sihat::getAudioExcerpt(u.soundFile, firstPeakSample, firstPeakSample + settings.modeWindow);
    std::vector<Synth::ModalComponent> modes = extractModalComponents(audioSlice, settings.overFrames);

    Synth::TModes tmodes;
    tmodes.modes = modes;
    tmodes.startInd = firstPeakSample;
    tmodes.length = settings.modeWindow;

    // Sihat::STransientResults results {tRange, tRise, peakVal, slicedEnv, tCentroid, tFlatness, tBands, tHarmonics, rms, settings.rmsHopSize, settings.hopSize, tHarmonics.hopSize, settings.nfft, settings.nfft / 2 + 1};
    Sihat::STransientResults results{tRange, settings.rmsHopSize, settings.hopSize, tHarmonics.hopSize, settings.nfft, settings.nfft / 2 + 1, tSpec.numFrames, tRise, peakVal, tEnv, tCentroid, tFlatness, tBands, tmodes, tHarmonics, rms};

    return results;
}

int TransientAnalysis::findFirstPeak(const std::vector<Sihat::Sample>& env, int firstSample)
{
    int peakIndex = 0;
    float maxRms = env[0].value;

    for (int i = 0; i < env.size(); i++)
    {
        if (env[i].value > maxRms)
        {
            maxRms = env[i].value;
            peakIndex = i;
        }
    }

    int finalSample = peakIndex * settings.rmsHopSize * 2;

    std::vector<float> aSegment = Sihat::getAudioExcerpt(u.soundFile, firstSample, finalSample);
    std::vector<float> aDiff = Sihat::getDifferential(aSegment); 

    float maxRawPeak = 0.0f;
    for (float val : aSegment)
    {
        float absVal = std::abs(val);
        if (absVal > maxRawPeak) maxRawPeak = absVal;
    }

    float threshold = maxRawPeak * settings.firstPeakThreshold; 

    for (size_t i = 1; i < aDiff.size() - 1; i++) 
    {
        float slopeIn = aDiff[i];
        float slopeOut = aDiff[i + 1];

        bool isExtremum = (slopeIn >= 0.0f && slopeOut < 0.0f) || 
                          (slopeIn <= 0.0f && slopeOut > 0.0f);

        if (isExtremum)
        {
            float peakAmp = std::abs(aSegment[i]);

            if (peakAmp >= threshold)
            {
                return firstSample + static_cast<int>(i);
            }
        }
    }

    auto maxIt = std::max_element(aSegment.begin(), aSegment.end(), 
        [](float a, float b) { return std::abs(a) < std::abs(b); });
        
    return firstSample + std::distance(aSegment.begin(), maxIt);
}

std::vector<Synth::ModalComponent> TransientAnalysis::extractModalComponents(const std::vector<float>& audio, int numFrames)
{
    int N = static_cast<int>(audio.size());
    
    // --- SAFETY GUARDS ---
    if (N < 32) return {}; 
    if (N > 4096) N = 4096; // Hard cap to prevent memory overflow at 96kHz

    int numExpectedModes = settings.numModes;
    int L = N / 3; 
    int numRows = N - L + 1;

    // --- MATRIX PENCIL METHOD ---
    Eigen::MatrixXf H(numRows, L);
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < L; ++j) H(i, j) = audio[i + j];
    }

    Eigen::BDCSVD<Eigen::MatrixXf> svd(H, Eigen::ComputeFullV);
    Eigen::MatrixXf V = svd.matrixV(); 

    int K = 2 * numExpectedModes; 
    if (K > L - 1) K = L - 1; // Prevent out-of-bounds column access
    if (K <= 0) return {};

    Eigen::MatrixXf Vk = V.leftCols(K);
    Eigen::MatrixXf V1 = Vk.topRows(L - 1);
    Eigen::MatrixXf V2 = Vk.bottomRows(L - 1);
    Eigen::MatrixXf Z = V1.colPivHouseholderQr().solve(V2);

    Eigen::EigenSolver<Eigen::MatrixXf> es(Z);
    auto complexPoles = es.eigenvalues();

    std::vector<Synth::ModalComponent> components;
    float T = 1.0f / u.sampleRate; // Assuming 'u' is accessible in this scope

    for (int i = 0; i < complexPoles.size(); ++i)
    {
        std::complex<float> zk = complexPoles[i];
        float mag = std::abs(zk);
        float phaseAngle = std::arg(zk);

        // Filter valid poles (stable, positive frequencies only)
        if (mag <= 1.0f && phaseAngle > 0.0f)
        {
            Synth::ModalComponent comp;
            comp.decay = -std::log(mag) / T;
            comp.freq = phaseAngle / (2.0f * Sihat::PI * T);
            components.push_back(comp);
        }
    }

    // --- LEAST SQUARES (AMPLITUDE & PHASE) ---
    int P = static_cast<int>(components.size());
    if (P == 0) return components;

    // Build a real-valued matrix W to solve for A*cos() and B*sin()
    // Size is N rows by 2P columns (2 columns per mode)
    Eigen::MatrixXf W(N, 2 * P);
    Eigen::VectorXf x(N);

    for (int n = 0; n < N; ++n)
    {
        x(n) = audio[n];
        float time = n * T;

        for (int m = 0; m < P; ++m)
        {
            float w = 2.0f * Sihat::PI * components[m].freq;
            float decayEnv = std::exp(-components[m].decay * time);

            // W0 handles the Cosine component, W1 handles the Sine component
            W(n, 2 * m)     = decayEnv * std::cos(w * time);
            W(n, 2 * m + 1) = decayEnv * -std::sin(w * time); 
        }
    }

    // Solve W * beta = x
    Eigen::VectorXf beta = W.colPivHouseholderQr().solve(x);

    // Reconstruct the physical amplitude and phase from the solved weights
    for (int m = 0; m < P; ++m)
    {
        float beta0 = beta(2 * m);     // Cosine weight
        float beta1 = beta(2 * m + 1); // Sine weight

        // R = sqrt(A^2 + B^2)
        components[m].amp = std::sqrt(beta0 * beta0 + beta1 * beta1);
        
        // Phase = atan2(B, A)
        components[m].phase = std::atan2(beta1, beta0); 
    }

    return components;
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
        if (seedPeak.freq < 250.0 ||
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
    int numFrames,
    bool trackFreq) // Added argument
{
    Sihat::TransientHarmonics tHarmonics;
    std::vector<float> floor;
    std::vector<Sihat::TransientOvertone> trajectories(targetPeaks.size());
    for (size_t i = 0; i < targetPeaks.size(); ++i) {
        trajectories[i].targetOvertone = targetPeaks[i];
        trajectories[i].envelope.reserve(numFrames); 
    }

    // 1. Decouple hop size and NFFT
    int totalRangeSamples = range.endSample - range.initSample;
    int frameStep = totalRangeSamples / numFrames;
    if (frameStep < 1) frameStep = 1;
    tHarmonics.hopSize = frameStep;
    
    // Hardcoded NFFT for frequency resolution as requested
    int nfft = settings.overNfft; 
    int sr = static_cast<int>(u.sampleRate);
    const int Nout = nfft / 2 + 1;

    float* input = (float*)fftwf_malloc(sizeof(float) * nfft);
    fftwf_complex* output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * Nout);
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(nfft, input, output, FFTW_MEASURE);

    std::vector<float> mags(Nout);

    // Center of the first frame starts at initSample
    int startOffset = range.initSample - (nfft / 2) > 0 ? range.initSample - (nfft / 2) : 0;
    tHarmonics.startSample = startOffset;

    const float trackingMinCrest = 1.5f;

    for (int i = 0; i < numFrames; ++i) 
    {
        int winStart = startOffset + (i * frameStep);
        int winEnd = winStart + nfft;

        // Zero-padding logic
        std::fill_n(input, nfft, 0.0f);
        int copyStart = std::max(winStart, 0); 
        int copyEnd = std::min(winEnd, (int)audio.size());

        if (copyStart < copyEnd) {
            int dstOffset = copyStart - winStart;
            std::copy(audio.begin() + copyStart, audio.begin() + copyEnd, input + dstOffset);
        }

        // Apply Hann Window
        for (int w = 0; w < nfft; ++w) {
            float winVal = 0.5f * (1.0f - std::cos(Sihat::TWO_PI * w / (nfft - 1)));
            input[w] *= winVal;
        }

        fftwf_execute(plan);

        for (int k = 0; k < Nout; ++k) {
            double re = output[k][0];
            double im = output[k][1];
            mags[k] = std::sqrt(re * re + im * im);
        }

        // State trackers for the current frame to use in the dynamic noise floor
        std::vector<int> currentFrameBins(targetPeaks.size(), -1);
        std::vector<double> currentFrameMaxMags(targetPeaks.size(), -1.0);
        std::vector<double> currentFrameFreqs(targetPeaks.size(), 0.0);

        // --- PASS 1: Frequency Tracking (Dynamic or Static) ---
        for (size_t h = 0; h < targetPeaks.size(); ++h) 
        {
            float initialTargetFreq = targetPeaks[h].freq;
            if (initialTargetFreq <= 20.0f) continue; 

            int bestBin = -1;
            double maxMag = -1.0;
            double finalFreq = initialTargetFreq;

            if (trackFreq) 
            {
                // Dynamic Tracking
                float searchFreq = initialTargetFreq;
                if (!trajectories[h].envelope.empty()) {
                    searchFreq = trajectories[h].envelope.back().freq;
                }

                int k_center = std::round(searchFreq * nfft / sr);
                int k_spread = 4; // Widened to allow following fast pitch bends
                
                int k_min = std::max(1, k_center - k_spread);
                int k_max = std::min(Nout - 2, k_center + k_spread);

                for (int k = k_min; k <= k_max; ++k) {
                    if (mags[k] > mags[k - 1] && mags[k] > mags[k + 1]) { // Local max check
                        if (mags[k] > maxMag) {
                            maxMag = mags[k];
                            bestBin = k;
                        }
                    }
                }

                // Fallback if no local peak is found
                if (bestBin == -1) {
                    bestBin = k_center;
                    maxMag = mags[bestBin];
                }

                // Parabolic Interpolation
                double m1 = mags[bestBin - 1], m0 = mags[bestBin], p1 = mags[bestBin + 1];
                double denom = (m1 - 2 * m0 + p1);
                double delta = 0.0;
                
                if (std::fabs(denom) >= 1e-12) {
                    delta = 0.5 * (m1 - p1) / denom;
                }
                
                finalFreq = (bestBin + delta) * sr / double(nfft);
            }
            else 
            {
                // Static Tracking
                bestBin = std::round(initialTargetFreq * nfft / sr);
                
                // Safety bound check
                bestBin = std::max(0, std::min(bestBin, Nout - 1)); 
                
                maxMag = mags[bestBin];
                finalFreq = initialTargetFreq;
            }

            currentFrameBins[h] = bestBin;
            currentFrameMaxMags[h] = maxMag;
            currentFrameFreqs[h] = finalFreq / f0; // Store as ratio to fundamental for later use in the envelope
        }

        // --- PASS 2: Dynamic Noise Floor Calculation ---
        float totalMagSum = 0.0;
        int binsCounted = 0;
        const int exclusionSpread = 2; // Exclude bestBin +/- 2 bins to handle Hann window leakage

        for (int bin = 0; bin < Nout; ++bin) {
            bool isHarmonicLeakage = false;
            for (int hBin : currentFrameBins) {
                if (hBin != -1 && std::abs(bin - hBin) <= exclusionSpread) {
                    isHarmonicLeakage = true;
                    break;
                }
            }

            if (!isHarmonicLeakage) {
                totalMagSum += mags[bin];
                binsCounted++;
            }
        }

        float averageFloor = 1e-12;
        if (binsCounted > 0) averageFloor = totalMagSum / binsCounted;
        if (averageFloor < 1e-12) averageFloor = 1e-12;

        // Window Gain Compensation (* 2.0f) for the noise floor
        float ampFloor = Sihat::mag_to_amp(averageFloor, nfft) * 2.0f;
        floor.push_back(ampFloor);

        // --- PASS 3: Build Final Envelopes ---
        for (size_t h = 0; h < targetPeaks.size(); ++h) 
        {
            Sihat::TrackedPoint point;
            float initialTargetFreq = targetPeaks[h].freq;

            if (initialTargetFreq <= 20.0f) {
                point.freq = initialTargetFreq;
                point.crestFactor = 0.0f;
                point.amp = 0.0f;
                point.active = false;
                trajectories[h].envelope.push_back(point);
                continue;
            }

            int bestBin = currentFrameBins[h];
            double maxMag = currentFrameMaxMags[h];

            point.freq = currentFrameFreqs[h];
            point.crestFactor = maxMag / averageFloor;
            
            // Window Gain Compensation (* 2.0f) for the overtone
            point.amp = Sihat::ampToDb(Sihat::mag_to_amp(maxMag, nfft) * 2.0f); 

            // We rely entirely on crest factor to gate the active state now.
            point.active = (point.crestFactor >= trackingMinCrest);

            // Safety check: If inactive, let the tracker rest at the last known good frequency
            // rather than wildly snapping to background noise bins
            if (!point.active) {
                point.freq = (!trajectories[h].envelope.empty()) ? trajectories[h].envelope.back().freq : initialTargetFreq;
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