#pragma once

#include<string>
#include"include/ResynthHeader.h"
#include"include/ResynthConfig.h"

class Resynthesizer
{
public:
    Resynthesizer(const Synth::Sihat& sh, const ResynthConfig& con) : sihat(sh), config(con) {};
    ~Resynthesizer(){};

    void resynthesize();

private:
    const Synth::Sihat& sihat;
    const ResynthConfig& config;

    inline std::vector<float> sumVectors(
    const std::vector<float>& base, 
    const std::vector<float>& addition, 
    const int baseStartSample, 
    const int addStartSample = 0,
    const float baseMultiplier = 1.0,
    const float addMultiplier = 1.0)
    {
        // 1. Guard against negative offsets
        if (baseStartSample < 0 || addStartSample < 0) {
            return {}; // Return empty vector or handle error as appropriate
        }

        // 2. Calculate the valid span of the addition vector
        int additionSamples = 0;
        if (addStartSample < static_cast<int>(addition.size())) {
            additionSamples = static_cast<int>(addition.size()) - addStartSample;
        }

        // 3. Determine the final size of the output vector
        // The base span ends at: base.size() + baseStartSample
        // The addition span ends at: additionSamples + baseStartSample
        int finalBaseEnd = static_cast<int>(base.size()) + baseStartSample;
        int finalAddEnd = additionSamples + baseStartSample;
        int finalSize = std::max(finalBaseEnd, finalAddEnd);

        // 4. Allocate the new vector, initialized to silence (0.0f)
        std::vector<float> result(finalSize, 0.0f);

        // 5. Copy the base vector into its offset position
        if (!base.empty()) {
            std::copy(base.begin(), base.end(), result.begin() + baseStartSample);
        }

        // 6. Accumulate (mix) the addition vector into its position
        for (int i = 0; i < additionSamples; ++i) {
            result[baseStartSample + i] *= baseMultiplier;
            result[baseStartSample + i] += addition[addStartSample + i] * addMultiplier;
        }

        return result;
    }




    std::vector<float> genTransient();
    std::vector<float> getTransientHarmonics();
    std::vector<float> getTransientNoise();
    std::vector<float> getTransientNoise(int hopSize, int nfft);
    std::vector<float> applyEnvelope(const std::vector<float>& vect, const Synth::Envelope& env, const int startSample = 0);
    std::vector<float> applyEnvelopeMatching(const std::vector<float>& vect, const Synth::Envelope& env, const int startSample = 0);
    std::vector<float> genHarmonics();
    std::vector<float> applyFlatness(const std::vector<float>& vect, const std::vector<float>& flatness, const int hopSize);
    void matchRms(std::vector<float>& vect, const float rms);
    void matchPeak(std::vector<float>& vect, const float peak, const int maxsize = 0);

};