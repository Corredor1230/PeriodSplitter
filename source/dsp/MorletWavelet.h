// #include "dsp/WTransform.h"

// class MorletWavelet : public WTransform
// {
// public:
//     using WTransform::WTransform;

//     std::vector<float> process(const std::vector<float>& frequencies) override{
//         std::copy(signal.begin(), signal.end(), inBuffer);
//         fftwf_execute(fwdPlan);

//         std::vector<float> scalogram;
//         scalogram.reserve(frequencies.size() * nfft);

//         for(float freq : frequencies)
//         {
//             std::vector<float> wave = wavelet(freq);

//             for (int k = 0; k < nfft; ++k)
//             {
//                 if (k < (nfft / 2 + 1))
//                 {
//                     float w = wave[k];
//                     invBuffer[k][0] = fDomainBuffer[k][0] * w; //real
//                     invBuffer[k][1] = fDomainBuffer[k][1] * w; //imaginaty
//                 }
//                 else
//                 {
//                     invBuffer[k][0] = 0.0f;
//                     invBuffer[k][1] = 0.0f;
//                 }
//             }

//             fftwf_execute(invPlan);

//             for (int i = 0; i < nfft; ++i)
//             {
//                 float re = invBuffer[i][0];
//                 float im = invBuffer[i][1];
//                 float mag = std::sqrt(re*re + im*im);

//                 scalogram.emplace_back(mag);
//             }

//             return scalogram;

//         }
//     }
// private:

//     const float omega0 = 6.0f;
//     std::vector<float> wavelet(float freq) override
//     {
//         std::vector<float> wave;
//         wave.reserve(nfft / 2 + 1);

//         float fc = omega0 / (Sitrano::PI * 2.0f);
//         float scale = (fc * fs) / freq;

//         for (int i = 0; i < (nfft / 2 + 1); ++i)
//         {
//             float omegak = (Sitrano::PI * 2.0f * i) / static_cast<float>(nfft);

//             float exp = scale * omegak - omega0;
//             float weight = 2.0f * std::exp(-0.5f * exp * exp);

//             wave.emplace_back(weight);
//         }
//         return wave;
//     }
// };

#include "dsp/WTransform.h"
#include <algorithm>
#include <vector>
#include <cmath>

class MorletWavelet : public WTransform
{
public:
    using WTransform::WTransform;

    // --- MAIN ENTRY POINT ---
    // Scans 'scanRes' scales, picks 'numPartials' winners, and extracts their envelopes.
    std::vector<Sitrano::VariableRatePartial> analyze(int numPartials, int scanRes, float minFreq, float maxFreq)
    {
        // 1. Run the Scanner (Fast Energy Calculation)
        // We find the best frequencies without doing the heavy iFFT yet.
        std::vector<float> bestFreqs = scanAndPickPeaks(scanRes, numPartials, minFreq, maxFreq);

        // 2. Extract Envelopes with Variable Rates
        std::vector<Sitrano::VariableRatePartial> results;
        results.reserve(bestFreqs.size());

        // Perform the standard Forward FFT once (if not already done)
        // Ensure inBuffer is loaded and fwdPlan executed
        std::copy(signal.begin(), signal.end(), inBuffer);
        fftwf_execute(fwdPlan);

        for (float freq : bestFreqs)
        {
            // A. Determine Hop Size (Downsampling Factor)
            // Rule of Thumb: Keep sample rate > Frequency / 4 for envelopes
            int hop = 1;
            if (freq < 100.0f)       hop = 256;
            else if (freq < 500.0f)  hop = 64;
            else if (freq < 2000.0f) hop = 16;
            else if (freq < 8000.0f) hop = 4;
            
            // B. Get the Full Resolution Envelope (The heavy lifting)
            std::vector<float> fullResEnvelope = extractEnvelope(freq);

            // C. Decimate (Downsample)
            Sitrano::VariableRatePartial partial;
            partial.frequency = freq;
            partial.hopSize = hop;
            partial.data.reserve(fullResEnvelope.size() / hop + 1);

            for (size_t t = 0; t < fullResEnvelope.size(); t += hop)
            {
                partial.data.push_back(fullResEnvelope[t]);
            }

            results.push_back(partial);
        }

        return results;
    }

    // Standard Process (Overload if you still need raw output)
    std::vector<float> process(const std::vector<float>& frequencies) override
    {
        std::copy(signal.begin(), signal.end(), inBuffer);
        fftwf_execute(fwdPlan);

        std::vector<float> scalogram;
        scalogram.reserve(frequencies.size() * nfft);

        for(float freq : frequencies)
        {
            std::vector<float> wave = wavelet(freq);

            for (int k = 0; k < nfft; ++k)
            {
                if (k < (nfft / 2 + 1))
                {
                    float w = wave[k];
                    invBuffer[k][0] = fDomainBuffer[k][0] * w; //real
                    invBuffer[k][1] = fDomainBuffer[k][1] * w; //imaginaty
                }
                else
                {
                    invBuffer[k][0] = 0.0f;
                    invBuffer[k][1] = 0.0f;
                }
            }

            fftwf_execute(invPlan);

            for (int i = 0; i < nfft; ++i)
            {
                float re = invBuffer[i][0];
                float im = invBuffer[i][1];
                float mag = std::sqrt(re*re + im*im);

                scalogram.emplace_back(mag);
            }

            return scalogram;

        }
    }

private:
    const float omega0 = 6.0f;

    // --- PHASE 1: THE SCANNER ---
    // Calculates total energy of 4000 bands using Frequency Domain only (Huge Optimization)
    std::vector<float> scanAndPickPeaks(int scanCount, int numToKeep, float minFreq, float maxFreq)
    {
        // 1. Generate Scan Frequencies (Logarithmic)
        std::vector<float> scanFreqs;
        scanFreqs.reserve(scanCount);
        float logMin = std::log10(minFreq);
        float logMax = std::log10(maxFreq);
        
        for (int i = 0; i < scanCount; ++i) {
            float t = (float)i / (float)(scanCount - 1);
            float f = std::pow(10.0f, logMin + (logMax - logMin) * t);
            scanFreqs.push_back(f);
        }

        // 2. Calculate Energy Profile
        // We use Parseval's Theorem: Sum of Energy in Freq Domain = Sum of Energy in Time Domain
        // This avoids running 4000 Inverse FFTs!
        std::vector<float> energyProfile(scanCount, 0.0f);

        // Ensure we have the signal spectrum ready
        std::copy(signal.begin(), signal.end(), inBuffer);
        fftwf_execute(fwdPlan);

        for (int i = 0; i < scanCount; ++i)
        {
            std::vector<float> wFilter = wavelet(scanFreqs[i]);
            float bandEnergy = 0.0f;

            // Multiply Signal * Wavelet and sum the squared magnitude
            for (int k = 0; k < (nfft / 2 + 1); ++k)
            {
                // Complex Signal: (a + ib)
                float rS = fDomainBuffer[k][0];
                float iS = fDomainBuffer[k][1];
                float w = wFilter[k];

                // Result = (rS*w) + i(iS*w)
                // Magnitude^2 = (rS*w)^2 + (iS*w)^2 = w^2 * (rS^2 + iS^2)
                float magSq = w * w * (rS * rS + iS * iS);
                bandEnergy += magSq;
            }
            energyProfile[i] = bandEnergy;
        }

        // 3. Peak Picking with Neighbor Suppression
        std::vector<float> winners;
        int suppressionWidth = scanCount / 50; // Suppress ~2% of the spectrum around a peak

        for (int k = 0; k < numToKeep; ++k)
        {
            auto maxIt = std::max_element(energyProfile.begin(), energyProfile.end());
            float maxVal = *maxIt;
            if (maxVal < 0.00001f) break; 

            int idx = std::distance(energyProfile.begin(), maxIt);
            winners.push_back(scanFreqs[idx]);

            // Nuke neighbors
            int start = std::max(0, idx - suppressionWidth);
            int end = std::min((int)energyProfile.size(), idx + suppressionWidth);
            for(int j = start; j < end; ++j) energyProfile[j] = 0.0f;
        }
        
        // Sort frequencies low-to-high for nicer output
        std::sort(winners.begin(), winners.end());
        return winners;
    }

    // --- PHASE 2: EXTRACTION ---
    // Performs the actual CWT (Forward -> Multiply -> Inverse -> Magnitude)
    std::vector<float> extractEnvelope(float freq)
    {
        std::vector<float> wave = wavelet(freq);
        std::vector<float> envelope;
        envelope.reserve(nfft);

        // Convolution in Frequency Domain
        for (int k = 0; k < nfft; ++k)
        {
            if (k < (nfft / 2 + 1)) {
                float w = wave[k];
                invBuffer[k][0] = fDomainBuffer[k][0] * w;
                invBuffer[k][1] = fDomainBuffer[k][1] * w;
            } else {
                invBuffer[k][0] = 0.0f;
                invBuffer[k][1] = 0.0f;
            }
        }

        fftwf_execute(invPlan);

        // Compute Magnitude
        for (int i = 0; i < nfft; ++i)
        {
            float re = invBuffer[i][0];
            float im = invBuffer[i][1];
            // Normalize by nfft if needed, but relative mag is fine usually
            envelope.push_back(std::sqrt(re*re + im*im));
        }
        return envelope;
    }

    // --- WAVELET GENERATION ---
    std::vector<float> wavelet(float freq) override
    {
        std::vector<float> wave;
        wave.reserve(nfft / 2 + 1);

        float fc = omega0 / (Sitrano::PI * 2.0f);
        float scale = (fc * fs) / freq;

        for (int i = 0; i < (nfft / 2 + 1); ++i)
        {
            float omegak = (Sitrano::PI * 2.0f * i) / static_cast<float>(nfft);
            float exponent = scale * omegak - omega0;
            // Note: Added factor 2.0 to compensate for analytic signal energy loss
            float weight = 2.0f * std::exp(-0.5f * exponent * exponent); 
            wave.emplace_back(weight);
        }
        return wave;
    }
};