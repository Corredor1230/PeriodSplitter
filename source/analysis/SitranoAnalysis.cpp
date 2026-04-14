#include"SitranoAnalysis.h"
#include"file/SihatFile.h"
#include"dsp/MedianFilter.h"

Sitrano::Results Analyzer::analyze(
    const Sitrano::AnalysisUnit& unit,
    const Sitrano::Settings& settings
)
{
    Sitrano::Results results;

    int provisionalHop = unit.sampleRate / 60;
    int provWindowNum = (unit.soundFile.size() - provisionalHop - mConfig.startSample) / provisionalHop;

    Sitrano::AnalysisUnit hUnit = unit;
    Sitrano::AnalysisUnit pUnit = unit;

    if (settings.sourceSeparation)
    {
        std::vector<float> harmonicAudio;
        std::vector<float> percussiveAudio;

        MedianFilter filter(mConfig.hpSettings);
        std::vector<float> interleavedHP = filter.processAudio(unit.soundFile);

        for (int i = 0; i < interleavedHP.size(); i++)
        {
            if (i % 2 == 0) harmonicAudio.push_back(interleavedHP[i]);
            else if (i % 2 == 1) percussiveAudio.push_back(interleavedHP[i]);
        }

        SihatFile::exportSeparatedAudio(interleavedHP, SihatFile::openFolderDialog(), (int)unit.sampleRate, "septest");

        hUnit.soundFile = harmonicAudio;
        pUnit.soundFile = percussiveAudio;
    }

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

    if (settings.transientSeparation)
    {
        Transient t(unit, mConfig.tSettings, mConfig.tfftSettings, results.pitch);
        results.tResults = t.findStartTransient();
    }
    else
    {
        results.tResults = { 0, 0 };
    }

    //Use the pitch to find the periods' zero crossings
    if (settings.periodAnalysis)
    {
        float pitch = results.pitch;
        if (!(pitch > mConfig.pSettings.minFreq)) pitch = unit.sampleRate / provisionalHop;

        PeriodCutter cutter{ unit, mConfig.cSettings, results.pitch, results.tResults.range };
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
        int firstSample = 0;

        if (!mConfig.oSettings.chooseFirstSample) firstSample = results.tResults.range.initSample;
        else firstSample = mConfig.oSettings.overtoneFirstSample;

        for (int i = 0; i < mConfig.oSettings.fftSize; ++i)
        {
            checkSignal.push_back(unit.soundFile[i + firstSample]);
        }
        OvertoneFinder finder{ unit, mConfig };
        results.topFreqs = finder.getRelevantOvertones(checkSignal, results.pitch);
        std::vector<float> freqs;
        for (int i = 0; i < results.topFreqs.size(); i++)
        {
            freqs.push_back(results.topFreqs[i].freq);
        }
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
        HarmonicTracker tracker(unit, mConfig, results.topFreqs, results.sampleList, results.tResults.range.endSample);
        results.hResults = tracker.analyze();
    }

    if (settings.noiseAnalysis)
    {
        NoiseTracker tracker(mConfig.nSettings, unit, results, 20.f);
        results.noise = tracker.analyze();
    }



    return results;
}