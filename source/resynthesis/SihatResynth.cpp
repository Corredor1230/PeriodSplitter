#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include"SihatResynth.h"
#include"file/SihatFile.h"

void Resynthesizer::resynthesize(const ResynthConfig& config)
{

    std::vector<float> harmonics;
    std::vector<float> transient;
    std::vector<float> output;

    //std::vector<float> transient = genTransient();
    if (config.resynthHarmonics)
    {harmonics = genHarmonics();}
    //std::vector<float> output = sumVectors(harmonics, transient, 0, 0);

    SihatFile::writeAudioFile(harmonics, config.filename + config.extension, config.r_outDir, sihat.header.sampleRate);
    
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