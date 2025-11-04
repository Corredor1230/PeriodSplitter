#include<sndfile.h>
#include<fftw3.h>
#include<iostream>
#include<math.h>
#include<string.h>
#include<vector>
#include<Windows.h>
#include<commdlg.h>
#include<filesystem>
#include"dsp/FIRFilter.h"
#include"gui/GUI.h"
#include"support/Splitter.h"
#include"analysis/SitranoAnalysis.h"
#include"include/csv.h"

namespace fs = std::filesystem;

std::string openFileDialog() {
    OPENFILENAMEA ofn;
    char szFile[260]; //Buffer for file name
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFile = szFile;
    ofn.lpstrFile[0] = '\0';
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "WAV Files\0*wav\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = nullptr;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = nullptr;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn) == TRUE) {
        return std::string(ofn.lpstrFile);
    }
    return "";
}

std::vector<std::string> getFileListFromExtension(const std::string& inDirectory, const std::string& extension)
{
    std::vector<std::string> filenames;

    std::string directory;

    if (inDirectory == "media")
    {
        directory = "C:/Users/usuario/Documents/Programming/CMake_Learning/PeriodicSplitter/media";
    }
    else if (inDirectory == "source")
    {
        directory = "C:/Users/usuario/Documents/Universidad/Tokyo_geijutsu_daigaku/2025_01/Master_Thesis/Media/Elec_Guitar/splitAudios";
    }
    else
    {
        directory = "C:/Users/usuario/Documents/Programming/CMake_Learning/PeriodicSplitter/media";
    }

    try {
        for (const auto& entry : fs::directory_iterator(directory)) {
            if (entry.is_regular_file() && entry.path().extension() == extension) {
                filenames.push_back(entry.path().filename().string());
            }
        }
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << '\n';
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }

    return filenames;
}

std::vector<float> getAudioFromFile(const std::string& filename, SF_INFO& sfInfo, float maxLength = 30.f)
{
    SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &sfInfo);
    int maxNumSamples = sfInfo.samplerate * maxLength;
    const int numFrames = sfInfo.frames;
    int numSamps = sfInfo.frames > maxNumSamples ? maxNumSamples : sfInfo.frames;
    std::vector<std::vector<float>> delaced(sfInfo.channels, std::vector<float>(numSamps));
    std::vector<float> samples(sfInfo.frames * sfInfo.channels);
    sf_readf_float(file, samples.data(), sfInfo.frames);

    for (int chan = 0; chan < sfInfo.channels; chan++)
    {
        for (int samp = 0; samp < numSamps; samp++)
        {
            delaced[chan][samp] = samples[samp * sfInfo.channels + chan];
        }
    }

    return delaced[0];
}

int main()
{

    bool bulkProcess = true;

    //General settings
    int maxHarmonics    = 32;
    int N               = 16384 * 2;
    int hopSize         = 1024;
    int startSample     = 0;
    float tolerance     = 100.0;
    float hTolerance    = 25.0;
    bool verbose        = true;
    Sitrano::WindowStyle w = Sitrano::WindowStyle::audioChunk;

    // Pitch config
    float modeThreshold = 3.0; //This is for finding the mode of the pitch array
    float tInCents      = 50.f; //Tolerance to compare found pitch with metadata
    float minFreq       = 60.f;
    float maxFreq       = 1300.f;

    //Transient config
    int tStartSample    = 0;
    bool tUseMs         = false;
    int rmsSampleSize   = 128;
    int rmsHop          = rmsSampleSize / 2;
    float transientRms  = 1.0f; //Measured in milliseconds
    float rmsHopRatio   = 1.0f;
    float tFactor       = 3.0f; //How big of an increase between RMS windows is a transient
    float tThreshold    = 0.1f;
    int preAttack       = 20;

    //Correlation config
    float periodOffset  = 50.0f; //How far after the transient should the analysis start
    float corrThreshold = 0.75f; //How similar should adjacent periods be

    //Overtone config
    bool overtoneTolerance  = true;
    bool oTolerance         = 200.f;
    bool postTransientStart = true;
    int overtoneFirstSample = 2000;
    bool useCustomSignal    = true;
    bool sumAmplitudes      = true;
    int oNfft               = N * 2;
    float ignoreThreshold   = -40.f;
    bool setAbsThreshold    = false;

    std::vector<std::string> fileList = getFileListFromExtension("source", ".wav");

    std::string filename = openFileDialog();

    SF_INFO sfInfo;
    std::vector<float> delaced = getAudioFromFile(filename, sfInfo);
    Sitrano::normalizeByMaxAbs(delaced);

    //CSV creation
    std::string csvName = Sitrano::getRawFilename(filename);

    Sitrano::AnalysisUnit unit{ 
        delaced, 
        filename,
        (float)sfInfo.samplerate
    };

    int overtoneStart = (int)((double)unit.soundFile.size() * 0.12);
    Sitrano::Settings settings{
        true,
        true,
        true,
        true,
        true,
        true
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
        Analyzer ana(config);

        //a.analyze(settings);
        Sitrano::Results r = ana.analyze(unit, settings);

        std::vector<float> singleFreqs;

        for (int i = 0; i < r.hResults.freqs.size(); i++) {
            if (r.hResults.freqs[i].size() <= 0) continue;
            singleFreqs.push_back(r.hResults.freqs[i][0]);
        }

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
        // --- Create the Analyzer ONCE, outside the loop ---
        Analyzer ana(config);

        // --- Get the list of files to process ---
        std::vector<std::string> fileList = getFileListFromExtension("source", ".wav");

        if (fileList.empty()) {
            std::cerr << "No .wav files found." << std::endl;
            return -1;
        }

        std::cout << "Found " << fileList.size() << " files to process." << std::endl;

        // --- This is your new batch processing loop ---
        for (const std::string& filename : fileList)
        {
            // Add a try-catch block for robust error handling
            try
            {
                std::cout << "--- Processing: " << filename << " ---" << std::endl;

                // --- All the per-file logic is moved inside the loop ---
                SF_INFO sfInfo;
                std::vector<float> delaced = getAudioFromFile(filename, sfInfo);
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

                // --- All your post-processing (filtering, etc.) ---
                std::vector<float> singleFreqs;
                for (int i = 0; i < r.hResults.freqs.size(); i++) {
                    if (r.hResults.freqs[i].size() <= 0) continue;
                    singleFreqs.push_back(r.hResults.freqs[i][0]);
                }

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
                    r.hResults.freqs, r.pitch, "sihat", "DATA_" + csvName, ".sihat");

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
