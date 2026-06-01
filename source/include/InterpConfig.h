#pragma once
#include "include/ResynthHeader.h"
#include "include/SitranoHeader.h"

namespace Interpolation
{

    enum class InterpType : uint32_t
    {
        Linear = 0,
        Spline = 1,
        Hybrid = 2
    };

    struct InterpParameters
    {
        bool hFreq = true;
        bool hAmp = true;
        bool tBands = true;
        bool riseTime = true;
        bool flatness = true;
    };

    struct InterpConfig
    {
        InterpType type = InterpType::Linear;
        float alpha = 0.5;
    };

    inline std::vector<float> f_interpolate(const std::vector<float>& vec1, const std::vector<float>& vec2)
    {
        int maxsize = std::max(vec1.size(), vec2.size());
        std::vector<float> out(maxsize, 0.0);

        
    }
}