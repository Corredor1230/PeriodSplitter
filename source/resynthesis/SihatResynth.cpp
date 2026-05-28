#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include"SihatResynth.h"
#include"file/SihatFile.h"

void Resynthesizer::resynthesize()
{

    std::vector<float> harmonics;
    std::vector<float> transient;
    std::vector<float> output;
    float harmAmp = Sihat::dbToAmp(config.harmOutDB);
    float transAmp = Sihat::dbToAmp(config.transOutDB);

    //std::vector<float> transient = genTransient();
    if (config.resynthHarmonics)
    {
        harmonics = genHarmonics();
        std::vector<float> temp(sihat.transient.meta.tEnd, 0.0);
        output = sumVectors(temp, harmonics, sihat.transient.riseTime * 2);
    }

    if (config.resynthTPercussive || config.resynthTHarmonic)
    {
        transient = genTransient();
        int maxsize = std::max(transient.size(), output.size());
        std::vector<float> temp(maxsize, 0.0);
        if (output.size() > 0) output = sumVectors(output, transient, 0, 0, harmAmp, transAmp);
        else output = sumVectors(temp, transient, 0);
    }

    matchPeak(output, 0.6, 96000);

    SihatFile::writeAudioFile(output, config.filename + config.extension, config.r_outDir, sihat.header.sampleRate);
    
}

std::vector<float> Resynthesizer::genHarmonics() {
    const auto& harm = sihat.harmonic;
    float sampleRate = static_cast<float>(sihat.header.sampleRate);
    size_t numOvertones = harm.amp.size();
    size_t numFrames = harm.indices.size();

    if (numFrames == 0 || numOvertones == 0) return {};

    // Define our absolute starting sample
    uint32_t startSample = sihat.transient.meta.tStart; 

    uint32_t lastHopSize = (numFrames > 1) ? 
        (harm.indices.back() - harm.indices[numFrames - 2]) : 1024;
    size_t totalSamples = harm.indices.back() + lastHopSize;

    // If our start point is past the end of our data, return silence
    if (startSample >= totalSamples) {
        return std::vector<float>(totalSamples, 0.0f);
    }

    std::vector<float> output(totalSamples, 0.0f);
    const float twoPi = 2.0f * Synth::PI;

    for (size_t k = 0; k < numOvertones; ++k) {
        // Phase resets to 0 for every overtone, guaranteeing they 
        // all start perfectly aligned at 'startSample'
        //if (k != 7) continue;
        float currentPhase = 0.0f;
        if (harm.amp[k].size() <= 0) continue;

        // Iterate through all frames except the last one
        for (size_t f = 0; f < numFrames - 1; ++f) {
            uint32_t frameStart = harm.indices[f];
            uint32_t frameEnd = harm.indices[f + 1];

            // Skip frames entirely if they happen before our start sample
            if (frameEnd <= startSample) continue;

            float startAmp = harm.amp[k][f];
            float endAmp = harm.amp[k][f + 1];
            float startFreq = harm.freq[k][f];
            float endFreq = harm.freq[k][f + 1];

            uint32_t currentWriteIndex = frameStart;

            // NEW: If this frame contains our startSample, we are starting mid-frame.
            if (frameStart < startSample) {
                currentWriteIndex = startSample;
                uint32_t originalFrameLength = frameEnd - frameStart;
                
                // Find exactly how far into the frame we are (0.0 to 1.0)
                float ratio = static_cast<float>(startSample - frameStart) / originalFrameLength;

                // Interpolate exact starting values so there are no sudden jumps
                startAmp = startAmp + ratio * (endAmp - startAmp);
                startFreq = startFreq + ratio * (endFreq - startFreq);
            }

            uint32_t synthesisLength = frameEnd - currentWriteIndex;
            if (synthesisLength == 0) continue;

            // Calculate deltas for the remainder of the frame
            float ampDelta = (endAmp - startAmp) / synthesisLength;
            float freqDelta = (endFreq - startFreq) / synthesisLength;

            float currentAmp = startAmp;
            float currentFreq = startFreq;

            for (uint32_t i = 0; i < synthesisLength; ++i) {
                output[currentWriteIndex + i] += currentAmp * std::sin(currentPhase);
                currentPhase += twoPi * (currentFreq / sampleRate);
                if (currentPhase >= twoPi) currentPhase -= twoPi;

                currentAmp += ampDelta;
                currentFreq += freqDelta;
            }
        }
        
        // Handle the final frame (fade out)
        uint32_t finalStart = harm.indices.back();
        uint32_t finalEnd = finalStart + lastHopSize;

        if (finalEnd > startSample) {
            uint32_t currentWriteIndex = finalStart;
            float currentAmp = harm.amp[k].back();
            float currentFreq = harm.freq[k].back();
            float endAmp = 0.0f; // Fade to absolute 0

            // Edge case: what if startSample is inside the fade-out tail?
            if (finalStart < startSample) {
                currentWriteIndex = startSample;
                float ratio = static_cast<float>(startSample - finalStart) / lastHopSize;
                currentAmp = currentAmp + ratio * (endAmp - currentAmp);
            }

            uint32_t synthesisLength = finalEnd - currentWriteIndex;
            float ampDelta = (endAmp - currentAmp) / synthesisLength;

            for (uint32_t i = 0; i < synthesisLength; ++i) {
                output[currentWriteIndex + i] += currentAmp * std::sin(currentPhase);
                currentPhase += twoPi * (currentFreq / sampleRate);
                if (currentPhase >= twoPi) currentPhase -= twoPi;

                currentAmp += ampDelta;
            }
        }
    }

    // 3. RMS Normalization (Only calculate over the active synthesized area)
    uint32_t activeSamples = sihat.header.sampleRate;
    double sumSquares = 0.0;
    
    // Start measuring only from where the audio actually begins
    for (uint32_t i = startSample; i < startSample + sihat.header.sampleRate; ++i) {
        sumSquares += output[i] * output[i];
    }
    
    float currentRms = std::sqrt(sumSquares / activeSamples);
    
    if (currentRms > 1e-8f && harm.rms > 0.0f) {
        float scalar = harm.rms / currentRms;
        for (uint32_t i = startSample; i < totalSamples; ++i) {
            output[i] *= scalar;
        }
    }

    return output;
}

std::vector<float> Resynthesizer::genTransient() {

    int tLen = sihat.transient.meta.tEnd;
    int first = sihat.transient.meta.tStart;

    float percDB = Sihat::dbToAmp(config.percDB);
    float harmDB = Sihat::dbToAmp(config.harmDB);

    std::vector<float> out(tLen, 0.0);

    if (config.resynthTHarmonic)
    {
        std::vector<float> harm = getTransientHarmonics();
        out = sumVectors(out, harm, first, 0, 1.0, harmDB);
    }

    if (config.resynthTPercussive)
    {
        std::vector<float> noise = getTransientNoise(sihat.transient.meta.specHopSize, sihat.transient.meta.specNfft);
        noise = applyEnvelopeMatching(noise, sihat.transient.envelope, 0);
        noise = applyEnvelope(noise, sihat.transient.envelope);
        //noise = applyFlatness(noise, sihat.transient.flatness, sihat.transient.meta.specHopSize);
        //matchRms(noise, sihat.transient.rms);
        out = sumVectors(out, noise, first, 0, 1.0, percDB);
    }

    out = applyEnvelope(out, sihat.transient.envelope, first);

    float outrms = Sihat::getRmsValue(out, 0, out.size());
    float ratio = sihat.transient.rms / outrms;
    for (auto& o : out) {
        o *= (ratio == 0.0) ? 1.0 : ratio;
    }

    return out;
}

std::vector<float> Resynthesizer::getTransientHarmonics() {
    const auto& transient = sihat.transient;

    const uint32_t sampleRate = sihat.header.sampleRate;
    const uint32_t hopSize = transient.meta.harmHopSize;
    const uint32_t numSamples = transient.meta.tEnd - transient.meta.tStart;

    std::vector<float> output(numSamples, 0.0);

    for (const auto& overtone : transient.overtones)
    {
        const auto& env = overtone.envelope;
        if (env.size() < 2) continue;

        double phase = overtone.target.phase;

        for (int frame = 0; frame < env.size() -1; frame++)
        {
            const auto& p0 = env[frame];
            const auto& p1 = env[frame + 1];

            const uint32_t startSample = frame * hopSize;
            const uint32_t endSample = (frame + 1) * hopSize;

            if (startSample >= numSamples) break;

            const uint32_t clampedEnd = std::min(endSample, numSamples);

            for (int n = startSample; n < clampedEnd; n++)
            {
                const float t = float(n - startSample) / float(hopSize);
                const float crest = p0.crestFactor + t * (p1.crestFactor - p0.crestFactor);
                const float freq = p0.freq + t * (p1.freq - p0.freq);
                const float amp = p0.amp + t * (p1.amp - p0.amp);

                output[n] += amp * std::sin(phase);

                phase += 2.0 * Synth::PI * freq / double(sampleRate);

                if (phase >= 2.0 * Synth::PI) phase -= 2.0 * Synth::PI;
            }
        }
    }

    //Pre-echo protection. This should improve performance significantly
    if (!transient.envelope.env.empty() && transient.riseTime > 0)
    {
        const auto& env = transient.envelope.env;
        const uint32_t envHop = transient.envelope.hopSize;
        
        // Find the peak of the original envelope to normalize it correctly
        float maxEnv = 0.0f;
        for (float val : env) { maxEnv = std::max(maxEnv, val); }
        if (maxEnv < 1e-6f) maxEnv = 1.0f;

        // We strictly shape only the attack phase to prevent double-decaying the tail
        const uint32_t attackSamples = std::min(static_cast<uint32_t>(transient.riseTime), numSamples);
        
        for (uint32_t n = 0; n < attackSamples; n++)
        {
            // Map the current time-domain sample 'n' to its fractional position in the envelope
            float envPos = static_cast<float>(n) / static_cast<float>(envHop);
            int frame = static_cast<int>(envPos);
            float t = envPos - frame;

            float e0 = (frame < env.size()) ? env[frame] / maxEnv : 1.0f;
            float e1 = (frame + 1 < env.size()) ? env[frame + 1] / maxEnv : 1.0f;

            // Linearly interpolate the exact target amplitude for this specific sample
            float targetAmp = e0 + t * (e1 - e0);

            // Clamp down the STFT pre-echo
            output[n] *= targetAmp;
        }
    }

    return output;
}

// Helper function for linear interpolation of the band envelope
float getInterpolatedBandAmp(const std::vector<float>& bands, int frame, int bin, int numBins, int numBands, int numFrames) {
    float binPosition = (float)bin / numBins * numBands; // Continuous position 0 to numBands
    int b0 = std::clamp(static_cast<int>(std::floor(binPosition - 0.5f)), 0, numBands - 1);
    int b1 = std::clamp(b0 + 1, 0, numBands - 1);
    
    float t = (binPosition - 0.5f) - b0;
    if (t < 0.0f) t = 0.0f; // Handle edge case at bin 0

    float energy0 = bands[b0 * numFrames + frame]; // Note: pass numFrames if needed in scope
    float energy1 = bands[b1 * numFrames + frame];
    
    // Interpolate energy, then convert to amplitude
    float interpEnergy = energy0 + t * (energy1 - energy0);
    return std::sqrt(interpEnergy / (numBins / numBands)); 
}

std::vector<float> Resynthesizer::getTransientNoise(int hopSize, int nfft)
{
    int numBins = sihat.transient.meta.specNumBins;
    int numFrames = sihat.transient.meta.specNumFrames;
    int numBands = sihat.transient.meta.numBands;

    int totalSamples = (numFrames * hopSize) + nfft;
    std::vector<float> output(totalSamples, 0.0f);
    
    // 1. Generate Continuous Time-Domain White Noise
    std::vector<float> continuousNoise(totalSamples);
    for (int i = 0; i < totalSamples; ++i) {
        continuousNoise[i] = ((rand() / (float)RAND_MAX) * 2.0f) - 1.0f;
    }

    // FFTW setup
    float* timeFrame = (float*)fftwf_malloc(sizeof(float) * nfft);
    fftwf_complex* complexFreq = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * numBins);
    
    // We need both Forward and Inverse plans now
    fftwf_plan fftPlan = fftwf_plan_dft_r2c_1d(nfft, timeFrame, complexFreq, FFTW_MEASURE);
    fftwf_plan ifftPlan = fftwf_plan_dft_c2r_1d(nfft, complexFreq, timeFrame, FFTW_MEASURE);

    std::vector<float> window = Sihat::getHannWindow(nfft);

    for (int i = 0; i < numFrames; ++i)
    {
        int readOffset = i * hopSize;

        // 2. Window a frame of the continuous noise
        for (int n = 0; n < nfft; ++n) {
            timeFrame[n] = continuousNoise[readOffset + n] * window[n];
        }

        // 3. To Frequency Domain
        fftwf_execute(fftPlan);

        // 4. Shape the spectrum using smoothly interpolated bands
        for (int j = 0; j < numBins; ++j)
        {
            // Calculate continuous interpolated amplitude
            float bandAmp = getInterpolatedBandAmp(sihat.transient.bands, i, j, numBins, numBands, numFrames);

            // Apply the amplitude envelope to the noise's natural phase/magnitude
            // We don't overwrite the phase; we scale the existing complex vector
            complexFreq[j][0] *= bandAmp; 
            complexFreq[j][1] *= bandAmp; 
        }

        // 5. Back to Time Domain
        fftwf_execute(ifftPlan);

        // 6. Overlap-Add with Synthesis Window
        for (int n = 0; n < nfft; ++n)
        {
            // FFTW unscaled c2r requires division by nfft
            float sample = (timeFrame[n] / nfft) * window[n];
            output[readOffset + n] += sample;
        }
    }

    fftwf_free(complexFreq);
    fftwf_free(timeFrame);
    fftwf_destroy_plan(fftPlan);
    fftwf_destroy_plan(ifftPlan);

    return output;
}

// std::vector<float> Resynthesizer::getTransientNoise(int hopSize, int nfft)
// {
//     int numBins = sihat.transient.meta.specNumBins;
//     int numFrames = sihat.transient.meta.specNumFrames;
//     int numBands = sihat.transient.meta.numBands;

//     // Total output size based on frames and hop size
//     int totalSamples = (numFrames * hopSize) + nfft;
//     std::vector<float> output(totalSamples, 0.0f);
    
//     int binsPerBand = numBins / numBands;

//     // FFTW setup for Inverse FFT
//     fftwf_complex* complexNoise = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * numBins);
//     float* timeFrame = (float*)fftwf_malloc(sizeof(float) * nfft);
//     fftwf_plan ifftPlan = fftwf_plan_dft_c2r_1d(nfft, complexNoise, timeFrame, FFTW_MEASURE);

//     for (int i = 0; i < numFrames; ++i)
//     {
//         // 1. Generate White Noise Spectrum
//         // In the frequency domain, white noise has random phase and flat magnitude.
//         for (int j = 0; j < numBins; ++j)
//         {
//             // Random phase between -PI and PI
//             float phase = ((rand() / (float)RAND_MAX) * Sihat::TWO_PI) - Sihat::PI;
            
//             // 2. Map bin 'j' to its corresponding band 'b'
//             int b = std::min(j / binsPerBand, numBands - 1);
            
//             // Extract energy and convert to amplitude
//             float bandEnergy = sihat.transient.bands[b * numFrames + i];
//             float bandAmp = std::sqrt(bandEnergy / binsPerBand);

//             // 3. Shape the noise (Polar to Cartesian conversion)
//             complexNoise[j][0] = bandAmp * std::cos(phase); // Real
//             complexNoise[j][1] = bandAmp * std::sin(phase); // Imaginary
//         }

//         // 4. IFFT back to time domain
//         fftwf_execute(ifftPlan);

//         // 5. Overlap-Add with Synthesis Window
//         int writeOffset = i * hopSize;
//         for (int n = 0; n < nfft; ++n)
//         {
//             // Standard Hann synthesis window
//             float winVal = 0.5f * (1.0f - std::cos(Sihat::TWO_PI * n / (nfft - 1)));
            
//             // Note: FFTW c2r transforms are unscaled, so we divide by nfft
//             float sample = (timeFrame[n] / nfft) * winVal;
            
//             output[writeOffset + n] += sample;
//         }
//     }

//     fftwf_free(complexNoise);
//     fftwf_free(timeFrame);
//     fftwf_destroy_plan(ifftPlan);

//     return output;
// }

std::vector<float> Resynthesizer::applyEnvelope(
    const std::vector<float>& vect,
    const Synth::Envelope& env, const int startSample)
{
    std::vector<float> out(vect.size(), 0.0);
    int hop = env.hopSize;
    int numFrames = std::min((vect.size() - startSample) / hop, env.env.size());
    float maxval = Sihat::getMaxVal(env.env);
    if (maxval < 1e-4) maxval = 1.0;

    for (int fr = 0; fr < numFrames - 1; fr++)
    {
        float e0 = env.env[fr] / maxval;
        float e1 = env.env[fr + 1] / maxval;

        float t = (e1 - e0) / static_cast<float>(hop);

        for (int s = 0; s < hop; s++)
        {
            out[fr * hop + s + startSample] = vect[fr * hop + s + startSample] * (e0 + (t * static_cast<float>(s)));
        }
    }

    return out;
}

std::vector<float> Resynthesizer::applyEnvelopeMatching(
    const std::vector<float>& vect,
    const Synth::Envelope& env, const int startSample)
{
    std::vector<float> out(vect.size(), 0.0);
    int hop = env.hopSize;
    int numFrames = std::min((vect.size() - startSample) / hop, env.env.size());
    
    // Safety floor to prevent division-by-zero or explosive gain
    const float epsilon = 1e-5f; 
    
    // We need to calculate the actual local amplitude (RMS) of the generated vector
    std::vector<float> currentRms(numFrames, 0.0f);
    for (int fr = 0; fr < numFrames; fr++)
    {
        float sumSq = 0.0f;
        for (int s = 0; s < hop; s++)
        {
            float val = vect[fr * hop + s + startSample];
            sumSq += val * val;
        }
        currentRms[fr] = std::sqrt(sumSq / hop);
    }

    // Now apply the Envelope Matching (Difference Ratio)
    for (int fr = 0; fr < numFrames - 1; fr++)
    {
        // 1. Get the target amplitudes (Note: NOT normalized by maxval. We want absolute levels)
        float target0 = env.env[fr];
        float target1 = env.env[fr + 1];

        // 2. Calculate the "Difference Ratio" (Gain factor) for this frame and the next
        float gain0 = target0 / (currentRms[fr] + epsilon);
        float gain1 = target1 / (currentRms[fr + 1] + epsilon);

        // Optional: Cap the maximum allowed gain to prevent artifacts if 'vect' is abnormally quiet
        gain0 = std::min(gain0, 10.0f); 
        gain1 = std::min(gain1, 10.0f);

        // 3. Calculate the linear interpolation step for the gain
        float gainStep = (gain1 - gain0) / static_cast<float>(hop);

        // 4. Apply the smoothed ratio
        for (int s = 0; s < hop; s++)
        {
            float currentGain = gain0 + (gainStep * static_cast<float>(s));
            out[fr * hop + s + startSample] = vect[fr * hop + s + startSample] * currentGain;
        }
    }

    return out;
}

std::vector<float> Resynthesizer::applyFlatness(const std::vector<float>& vec, const std::vector<float>& flat, const int hopSize)
{
    int numFrames = flat.size();
    std::vector<float> out(vec.size(), 0.0);
    float maxval = Sihat::getMaxVal(flat);
    if (maxval < 1e-4) maxval = 1.0;

    for (int f = 1; f < numFrames - 1; f++)
    {
        float flat0 = flat[f] / maxval;
        float flat1 = flat[f + 1] / maxval;

        float t = (flat1 - flat0) / static_cast<float>(hopSize);

        for (int s = 0; s < hopSize; s++)
        {
            int curr = f * hopSize + s;
            out[curr] = vec[curr] * (flat0 + (t * s));
        }
    }

    return out;
}

void Resynthesizer::matchRms(std::vector<float>& vec, const float rms)
{
    float inRms = Sihat::getRmsValue(vec, 0, vec.size());
    float mRms = rms;
    if (inRms == 0.0) return;
    float ratio = mRms / inRms;

    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] *= ratio;
    }
}

void Resynthesizer::matchPeak(std::vector<float>& vec, const float peak, const int maxsize)
{
    if (vec.size() == 0) return;

    int max = vec.size();
    float maxval = 0.0;
    if (maxsize != 0) max = maxsize;
    
    for (int i = 0; i < max; i++)
    {
        if (vec[i] > maxval) maxval = vec[i];
    }
    
    float ratio = maxval >= 1.0 ? peak / maxval : maxval / peak;

    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] *= ratio;
    }
}