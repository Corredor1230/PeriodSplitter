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
#include"analysis/HarmonicTracker.h"
#include"analysis/NoiseTracker.h"
#include"include/SitranoHeader.h"

class Analyzer {
    Sitrano::AnalysisConfig mConfig;

public:

    /*Analyzer(const Sitrano::AnalysisUnit&,
        Sitrano::Results& r,
        bool applyHann = true,
        float toleranceVal = 300.0);*/

    Analyzer(const Sitrano::AnalysisConfig& config);

    ~Analyzer() {
        //fftwf_destroy_plan(plan);
        //fftwf_free(input);
        //fftwf_free(output);
    }

    //void analyze(const Sitrano::Settings& settings);

    Sitrano::Results analyze(
        const Sitrano::AnalysisUnit& unit,
        const Sitrano::Settings& settings);

private:

    /*std::vector<Sitrano::Peak> findTopSorted(
        const std::vector<double>& input,
        int N, double fs, int startPos,
        bool useCentsTolerance = true,
        double tolValue = 15.0,
        bool sumAmplitudesInCluster = true);*/

    //std::vector<float> window;

    //const Sitrano::AnalysisUnit& ana;
    //Sitrano::Results& results;

    //fftwf_plan plan;
    //float* input;
    std::vector<double> checker;
    //fftwf_complex* output;

    //void initFFTW();
    //void applyHann(double* data, int size);
};
