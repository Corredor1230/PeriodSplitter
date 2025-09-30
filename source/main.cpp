#include<sndfile.h>
#include<fftw3.h>
#include<iostream>
#include<math.h>
#include<cmath>
#include<algorithm>
#include<string.h>
#include<vector>
#include<Windows.h>
#include<commdlg.h>
#include<filesystem>
#include<libpyincpp.h>
#include<ranges>
#include"PeriodCutter.h"
#include"FIRFilter.h"
#include"GUI.h"
#include"Splitter.h"
#include"HarmonicTracker.h"
#include"SitranoAnalysis.h"
#include"csv.h"

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

float findMode(const std::vector<float>& data, float threshold = 1.0f,
    float minFreq = 60.0, float maxFreq = 1200.0)
{
    if (data.empty()) return NAN;

    // Step 1: Sort the data
    std::vector<float> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Step 2: Slide a window and track the group with the highest count
    int maxCount = 0;
    float modeSum = 0.0f;
    int modeCount = 0;

    for (size_t i = 0; i < sortedData.size(); ++i) {
        float ref = sortedData[i];
        float sum = 0.0f;
        int count = 0;

        if (ref < minFreq || ref > maxFreq)
            continue;

        // Include all values within [ref, ref + threshold)
        for (size_t j = i; j < sortedData.size() && sortedData[j] - ref < threshold; ++j) {
            sum += sortedData[j];
            ++count;
        }

        if (count > maxCount) {
            maxCount = count;
            modeSum = sum;
            modeCount = count;
        }
    }

    float mode = modeSum / modeCount;

    return modeCount > 0 ? mode : NAN;
}

float getPitchFromFilename(std::string& filename)
{
    std::string name = filename;
    std::vector<std::string> parts;
    size_t start = 0, pos;

    while ((pos = name.find('_', start)) != std::string::npos) {
        parts.push_back(name.substr(start, pos - start));
        start = pos + 1;
    }
    parts.push_back(name.substr(start));

    float outPitch = std::stof(parts[2]);

    return outPitch;
}

float findPitch(const std::vector<float>& signal, float sampleRate, std::string& filename)
{
    PyinCpp pitchDetector(sampleRate);

    std::vector<float> pitches = pitchDetector.feed(signal);

    float foundPitch = findMode(pitches);

    //std::cout << "Current pitch: " << pitch << "Hz\n";

    float metaPitch = getPitchFromFilename(filename);
    float pitch{ 0.0 };
    float tolerance = 3.0;

    if (foundPitch > metaPitch - tolerance && foundPitch < metaPitch + tolerance)
        pitch = foundPitch;
    else
        pitch = metaPitch;

    std::cout << "Current pitch: " << pitch << 'Hz\n';

    return pitch;
}

std::vector<float> getAudioFromFile(std::string& filename, SF_INFO& sfInfo, float maxLength = 30.f)
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
    int maxHarmonics = 32;
    int N = 16384 * 2;
    float tolerance = 200.f;

    std::vector<std::string> fileList = getFileListFromExtension("source", ".wav");

    std::string filename = openFileDialog();

    SF_INFO sfInfo;
    std::vector<float> delaced = getAudioFromFile(filename, sfInfo);
    Sitrano::normalizeByMaxAbs(delaced);
    int startSample = 0;

    //CSV creation
    //std::string csvName = Sitrano::getRawFilename(filename);
    //Splitter theSplitter(sfInfo);
    //theSplitter.writeCsvFile(periodStart, csvName + ".csv");

    /*auto windows = theSplitter.loadCSV(csvName + ".csv", delaced.size());
    WaveformViewer viewer(delaced, windows);
    viewer.run();*/

    Sitrano::Results r;
    Sitrano::AnalysisUnit unit{ 
        delaced, 
        filename,
        (float)sfInfo.samplerate, 
        maxHarmonics, 
        N,
        startSample};

    Sitrano::Settings settings{
        true,
        true,
        true,
        true,
        true
    };

    Analyzer a(
        unit, 
        r, 
        true, 
        300.0);

    a.analyze(settings);

    std::vector<float> singleFreqs;

    for (int i = 0; i < r.freqs.size(); i++) {
        singleFreqs.push_back(r.freqs[i][0]);
    }

    for (int i = 0; i < r.amps.size(); i++)
    {
        Sitrano::filterVector(r.amps[i], 4);
    }

    //theSplitter.writeCsvFile(r.amps, "AMP" + csvName + ".csv", singleFreqs, periodStart);
    //theSplitter.writeCsvFile(phases, "PHA" + csvName + ".csv", singleFreqs);

	return 0;
}
