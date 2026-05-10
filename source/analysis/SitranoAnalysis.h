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
#include"analysis/TransientAnalysis.h"
#include"include/SitranoHeader.h"

class Analyzer {
    const Sihat::AnalysisConfig mConfig;

public:

    Analyzer(const Sihat::AnalysisConfig& config) : mConfig(config) 
    {};

    ~Analyzer() {}

    Sihat::Results analyze(
        const Sihat::AnalysisUnit& unit,
        const Sihat::Settings& settings);

private:
    std::vector<double> checker;
};
