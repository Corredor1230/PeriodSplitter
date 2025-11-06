#include<sndfile.h>
#include<iostream>
#include<math.h>
#include<string.h>
#include<vector>
#include"dsp/FIRFilter.h"
#include"gui/GUI.h"
#include"support/Splitter.h"
#include"analysis/SitranoAnalysis.h"
#include"include/csv.h"
#include"file/SihatFile.h"

void processFile(const std::string& filename,
    const Sitrano::AnalysisConfig& config,
    const Sitrano::Settings& settings,
    const std::string& outDir = "")
{
    SF_INFO sfInfo;
    std::vector<float> delaced = SihatFile::getAudioFromFile(filename, sfInfo);
    Sitrano::normalizeByMaxAbs(delaced);

    Sitrano::AnalysisUnit unit{ delaced, filename, (float)sfInfo.samplerate };
    Analyzer ana(config);
    Sitrano::Results r = ana.analyze(unit, settings);

    for (auto& amps : r.hResults.amps)
        if (!amps.empty()) Sitrano::filterVector(amps, 4, true);
    for (auto& freqs : r.hResults.freqs)
        if (!freqs.empty()) Sitrano::filterVector(freqs, 40, false);

    std::string csvName = Sitrano::getRawFilename(filename);
    std::string dir = outDir.empty() ? "sihat" : outDir;

    Sitrano::saveHarmonicDataSihat(
        r.hResults.finalSamples, r.hResults.amps, r.hResults.freqs,
        r.pitch, dir, "DATA_" + csvName, ".sihat");
}

int main()
{
    bool bulkProcess = true;

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

    if (!bulkProcess)
    {
        std::string filename = SihatFile::openFileDialog();

        SF_INFO sfInfo;
        std::vector<float> delaced = SihatFile::getAudioFromFile(filename, sfInfo);
        Sitrano::normalizeByMaxAbs(delaced);

        //CSV creation
        std::string csvName = Sitrano::getRawFilename(filename);

        Sitrano::AnalysisUnit unit{
            delaced,
            filename,
            (float)sfInfo.samplerate
        };

        Analyzer ana(config);

        //a.analyze(settings);
        Sitrano::Results r = ana.analyze(unit, settings);

        for (int i = 0; i < r.hResults.amps.size(); i++)
        {
            if (r.hResults.amps[i].size() <= 0) continue;
            Sitrano::filterVector(r.hResults.amps[i], 4, true);
        }

        for (int i = 0; i < r.hResults.freqs.size(); i++) {
            if (r.hResults.freqs[i].size() <= 0) continue;
            Sitrano::filterVector(r.hResults.freqs[i], 40, false);
        }

        Sitrano::saveHarmonicDataSihat(r.hResults.finalSamples, r.hResults.amps,
            r.hResults.freqs, r.pitch, "sihat", "DATA_" + csvName, ".sihat");
    }

    else
    {
        // --- Get the list of files to process ---
        std::string bulkDir = SihatFile::openFolderDialog();
        std::string saveDir = SihatFile::openFolderDialog();
        SihatFile::DirAndFiles dnf = SihatFile::getFileListFromExtension(bulkDir, ".wav");

        if (dnf.files.empty()) {
            std::cerr << "No .wav files found." << std::endl;
            return -1;
        }

        std::cout << "Found " << dnf.files.size() << " files to process." << std::endl;

        // --- This is your new batch processing loop ---
        for (const std::string& fn : dnf.files)
        {
            SF_INFO sfInfo;
            std::string filename = dnf.directory + "/" + fn;
            std::vector<float> delaced = SihatFile::getAudioFromFile(filename, sfInfo);
            Sitrano::normalizeByMaxAbs(delaced);

            //CSV creation
            std::string csvName = Sitrano::getRawFilename(filename);

            Sitrano::AnalysisUnit unit{
                delaced,
                filename,
                (float)sfInfo.samplerate
            };

            Analyzer ana(config);

            // Add a try-catch block for robust error handling
            try
            {
                std::cout << "--- Processing: " << filename << " ---" << std::endl;

                // --- All the per-file logic is moved inside the loop ---
                SF_INFO sfInfo;
                std::vector<float> delaced = SihatFile::getAudioFromFile(filename, sfInfo);
                Sitrano::normalizeByMaxAbs(delaced);

                //CSV creation
                std::string csvName = Sitrano::getRawFilename(filename);

                Sitrano::AnalysisUnit unit{
                    delaced,
                    filename,
                    (float)sfInfo.samplerate
                };

                // Run analysis
                Sitrano::Results r = ana.analyze(unit, settings);

                for (int i = 0; i < r.hResults.amps.size(); i++)
                {
                    if (r.hResults.amps[i].size() <= 0) continue;
                    Sitrano::filterVector(r.hResults.amps[i], 4, true);
                }

                for (int i = 0; i < r.hResults.freqs.size(); i++) {
                    if (r.hResults.freqs[i].size() <= 0) continue;
                    Sitrano::filterVector(r.hResults.freqs[i], 40, false);
                }

                // --- Save the output file ---
                Sitrano::saveHarmonicDataSihat(r.hResults.finalSamples, r.hResults.amps,
                    r.hResults.freqs, r.pitch, saveDir, "DATA_" + csvName, ".sihat");

                std::cout << "--- Finished: " << filename << " ---" << std::endl;
            }
            catch (const std::exception& e)
            {
                // If one file fails, log it and continue to the next one
                std::cerr << "!!! FAILED to process " << filename << ": " << e.what() << std::endl;
            }
            catch (...)
            {
                std::cerr << "!!! FAILED to process " << filename << ": Unknown error" << std::endl;
            }
        }

        std::cout << "Batch processing complete." << std::endl;
    }

    return 0;
}
