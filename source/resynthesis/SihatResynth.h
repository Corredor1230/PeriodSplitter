#pragma once

#include<string>
#include"include/ResynthHeader.h"
#include"include/ResynthConfig.h"

class Resynthesizer
{
public:
    Resynthesizer(const Synth::Sihat& sh) : sihat(sh) {};
    ~Resynthesizer(){};

    void resynthesize(const ResynthConfig& config);

private:
    const Synth::Sihat& sihat;

    inline std::vector<float> sumVectors(
    const std::vector<float>& base, 
    const std::vector<float>& addition, 
    const int baseStartSample, 
    const int addStartSample = 0)
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
            result[baseStartSample + i] += addition[addStartSample + i];
        }

        return result;
    }




    std::vector<float> genTransient();
    std::vector<float> genHarmonics();

};