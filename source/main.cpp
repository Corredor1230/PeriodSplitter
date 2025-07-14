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
#include"Correlator.h"
#include"FIRFilter.h"
#include"GUI.h"
#include"Splitter.h"
#include"csv.h"

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

std::string getRawFilename(std::string& filename)
{
    std::string rawFilename;
    std::filesystem::path path = filename;
    rawFilename = path.filename().replace_extension("").string();

    return rawFilename;
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

float findPitch(const std::vector<float>& signal, float sampleRate)
{
    PyinCpp pitchDetector(sampleRate);

    std::vector<float> pitches = pitchDetector.feed(signal);

    float pitch = findMode(pitches);

    std::cout << "Current pitch: " << pitch << "Hz\n";

    return pitch;
}

void normalizeByMaxAbs(std::vector<float>& vec) 
{
    if (vec.empty()) return;

    // Find the maximum absolute value
    float maxAbs = 0.0f;
    for (float x : vec) {
        maxAbs = std::max<float>(maxAbs, std::abs(x));
    }

    if (maxAbs == 0.0f) {
        // Avoid division by zero — all elements are zero
        return;
    }

    // Normalize
    for (float& x : vec) {
        x /= maxAbs;
    }
}

int main()
{
    std::string filename = openFileDialog();

    SF_INFO sfInfo;
    SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &sfInfo);
    int maxNumSamples = sfInfo.samplerate * 10;
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

    int startSample = 0;
    float pitch = findPitch(delaced[0], sfInfo.samplerate);
    bool needsFiler = false;
    if (pitch > 250.0)
    {
        float cutoff = pitch > 300.0 ? pitch * 0.50 : 150.0;
        FIRHighPass filter(sfInfo.samplerate, cutoff, cutoff * 0.8);
        filter.processAudioFile(delaced[0], 1024);
    }
    normalizeByMaxAbs(delaced[0]);

    Correlator correlator(delaced[0], sfInfo);
    correlator.initialize(pitch);
    std::string csvName = getRawFilename(filename);

    std::vector<int> periodStart = correlator.getCorrelationZeroes();

    Splitter theSplitter(sfInfo);
    theSplitter.writeCsvFile(periodStart, csvName + ".csv");

    auto windows = theSplitter.loadCSV(csvName + ".csv", delaced[0].size());

    WaveformViewer viewer(delaced[0], windows);
    viewer.run();

	return 0;
}