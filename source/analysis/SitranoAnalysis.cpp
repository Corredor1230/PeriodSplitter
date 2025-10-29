#include"SitranoAnalysis.h"

Sitrano::Results Analyzer::analyze(
    const Sitrano::AnalysisUnit& unit,
    const Sitrano::Settings& settings
)
{
    Sitrano::Results results;

    int provisionalHop = unit.sampleRate / 60;
    int provWindowNum = (unit.soundFile.size() - provisionalHop - mConfig.startSample) / provisionalHop;

    //Find the pitch
    if (settings.pitchAnalysis)
    {
        PitchFinder finder{ unit, mConfig };
        results.pitch = finder.findPitch();
    }
    else
    {
        results.pitch = Sitrano::getPitchFromFilename(unit.filename);
    }

    //Use the pitch to find the periods' zero crossings
    if (settings.periodAnalysis)
    {
        float pitch = results.pitch;
        if (!pitch > mConfig.pConfig.minFreq) pitch = unit.sampleRate / provisionalHop;

        PeriodCutter cutter{ unit, mConfig.cSettings, results.pitch, mConfig.startSample };
        results.sampleList = cutter.findPeriodSamples();
    }
    else
    {
        for (int i = 0; i < provWindowNum; i++)
        {
            results.sampleList.push_back(i * provisionalHop + mConfig.startSample);
        }
    }

    //Check the most representative overtones
    if (settings.overtoneAnalysis)
    {
        std::vector<double> checkSignal;
        for (int i = 0; i < mConfig.oConfig.fftSize; ++i)
        {
            checkSignal.push_back(unit.soundFile[i + mConfig.oConfig.overtoneFirstSample]);
        }
        OvertoneFinder finder{ unit, mConfig };
        results.topFreqs = finder.getRelevantOvertones(checkSignal, results.pitch);
    }
    else
    {
        for (int i = 0; i < mConfig.numHarmonics; ++i)
        {
            Sitrano::Peak p{
                i * results.pitch,
                0.f,
                1.0 / (float)(i + 1)
            };
            results.topFreqs.push_back(p);
        }
    }


    //Check the amplitude envelope for each harmonic
    if (settings.harmonicAnalysis)
    {
        HarmonicTracker tracker(unit, mConfig, results.topFreqs, results.sampleList);
        results.hResults = tracker.analyze();
    }

    return results;
}