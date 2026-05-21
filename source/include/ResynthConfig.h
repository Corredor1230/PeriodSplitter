#pragma once

#include<string>
#include"include/SitranoHeader.h"

struct ResynthConfig{
    std::string filename = "gen";
    std::string extension = ".wav";
    std::string r_outDir = "./";
    std::string fullPath = r_outDir + filename + extension;
    bool resynthTransient = true;
    bool resynthHarmonics = true;
    bool separateOuts = false;
    bool changeOutLevel = false;
    float outLevelDB = 0.0;
    float outLevel = Sihat::dbToAmp(outLevelDB);
};