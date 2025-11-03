#pragma once

#include<vector>
#include<fftw3.h>
#include<cmath>
#include<math.h>
#include<algorithm>
#include<iostream>
#include"dsp/PitchFinder.h"
#include"dsp/PeriodCutter.h"
#include"analysis/OvertoneFinder.h"
#include"analysis/TransientRecognition.h"
#include"analysis/HarmonicTracker.h"
#include"analysis/NoiseTracker.h"
#include"include/SitranoHeader.h"

class Analyzer {
    const Sitrano::AnalysisConfig mConfig;

public:

    Analyzer(const Sitrano::AnalysisConfig& config) : mConfig(config) 
    {};

    ~Analyzer() {}

    Sitrano::Results analyze(
        const Sitrano::AnalysisUnit& unit,
        const Sitrano::Settings& settings);

private:
    std::vector<double> checker;
};
