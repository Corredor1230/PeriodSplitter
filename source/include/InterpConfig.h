#pragma once
#include "include/ResynthHeader.h"
#include "include/SitranoHeader.h"
#include "include/InterpHeader.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

namespace InterpConfig
{

    struct Parameters
    {
        bool hFreq = true;
        bool hAmp = true;
        bool tBands = true;
        bool riseTime = true;
        bool flatness = true;
    };

    struct Config
    {
        SihatInterpolation::Type type = SihatInterpolation::Type::Linear;
        std::string i_filename = "interpolation";
        uint32_t sampleRate = 96000;

        std::string pathA = "";
        std::string pathB = "";
        float alpha = 0.5;
        float f0 = 440.0;
    };
}