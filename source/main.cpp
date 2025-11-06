#include<iostream>
#include<string.h>
#include<vector>
#include"file/SihatFile.h"

int main()
{
    bool bulkProcess = true;

    std::string outDir      = "sihat";
    std::string prefix      = "DATA_";
    std::string extension   = ".sihat";

    //Analysis
    bool pitchAnalysis = true;
    bool transientSeparation = true;
    bool periodAnalysis = true;
    bool overtoneAnalysis = true;
    bool harmonicAnalysis = true;
    bool noiseAnalysis = true;

    //General settings
    int maxHarmonics = 32;
    int N = 16384 * 2;
    int hopSize = 1024;
    int startSample = 0;
    float tolerance = 100.0;
    float hTolerance = 25.0;
    bool verbose = true;
    Sitrano::WindowStyle w = Sitrano::WindowStyle::audioChunk;

    // Pitch config
    float modeThreshold = 3.0; //This is for finding the mode of the pitch array
    float tInCents = 50.f; //Tolerance to compare found pitch with metadata
    float minFreq = 60.f;
    float maxFreq = 1300.f;

    //Transient config
    int tStartSample = 0;
    bool tUseMs = false;
    int rmsSampleSize = 128;
    int rmsHop = rmsSampleSize / 2;
    float transientRms = 1.0f; //Measured in milliseconds
    float rmsHopRatio = 1.0f;
    float tFactor = 3.0f; //How big of an increase between RMS windows is a transient
    float tThreshold = 0.1f;
    int preAttack = 20;

    //Correlation config
    float periodOffset = 50.0f; //How far after the transient should the analysis start
    float corrThreshold = 0.75f; //How similar should adjacent periods be

    //Overtone config
    bool overtoneTolerance = true;
    bool oTolerance = 200.f;
    bool postTransientStart = true;
    int overtoneFirstSample = 2000;
    bool useCustomSignal = true;
    bool sumAmplitudes = true;
    int oNfft = N * 2;
    float ignoreThreshold = -40.f;
    bool setAbsThreshold = false;

    Sitrano::Settings settings{
        pitchAnalysis,
        transientSeparation,
        periodAnalysis,
        overtoneAnalysis,
        harmonicAnalysis,
        noiseAnalysis
    };

    Sitrano::OvertoneSettings oSettings{
        overtoneTolerance,
        oTolerance,
        postTransientStart,
        overtoneFirstSample,
        useCustomSignal,
        sumAmplitudes,
        oNfft,
        ignoreThreshold,
        setAbsThreshold
    };

    Sitrano::CorrelationSettings cSettings{
        periodOffset,
        corrThreshold
    };

    Sitrano::TransientSettings tSettings{
        tStartSample,
        tUseMs,
        rmsSampleSize,
        rmsHop,
        transientRms,
        rmsHopRatio,
        tFactor,
        tThreshold,
        preAttack
    };

    Sitrano::HarmonicSettings hSettings{
        true,
        w,
        hTolerance
    };

    Sitrano::PitchSettings pSettings{
        modeThreshold,
        tInCents,
        minFreq,
        maxFreq
    };

    Sitrano::AnalysisConfig config{
        maxHarmonics,
        N,
        hopSize,
        startSample,
        tolerance,
        pSettings,
        tSettings,
        cSettings,
        oSettings,
        hSettings,
        verbose
    };

    SihatFile::OutInfo info{
        outDir,
        prefix,
        extension
    };

    if (!bulkProcess)
    {
        std::string filename = SihatFile::openFileDialog();
        SihatFile::processFile(filename, info, config, settings);
    }
    else
    {
        std::string inputDir = SihatFile::openFolderDialog();
        std::string saveDir = SihatFile::openFolderDialog();
        SihatFile::processFolder(inputDir, info, config, settings);
    }

    return 0;
}
