#include<sndfile.h>
#include<fftw3.h>
#include<iostream>
#include<math.h>
#include<string.h>
#include<vector>
#include<Windows.h>
#include<commdlg.h>
#include"Correlator.h"
#include"Splitter.h"

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

int main()
{
    std::string filename = openFileDialog();

    SF_INFO sfInfo;
    SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &sfInfo);
    const int numFrames = sfInfo.frames;
    std::vector<std::vector<float>> delaced(sfInfo.channels, std::vector<float>(sfInfo.frames));
    std::vector<float> samples(sfInfo.frames * sfInfo.channels);
    sf_readf_float(file, samples.data(), sfInfo.frames);

    for (int chan = 0; chan < sfInfo.channels; chan++)
    {
        for (int samp = 0; samp < sfInfo.frames; samp++)
        {
            delaced[chan][samp] = samples[samp * sfInfo.channels + chan];
        }
    }

    int startSample = 0;

    Correlator correlator(delaced[0], sfInfo);

    std::string instrument("guit");

    correlator.calculateCorrelation(startSample);
    std::string corrFilename(instrument + "Correlation");
    //correlator.printCorrelation(corrFilename);
    std::string corrPeak(instrument + "Peak");
    //correlator.printCorrelationPeak(corrPeak);
    std::string zeroFilename(instrument + "Zeroes");
    //correlator.printZeroList(zeroFilename);

    Splitter theSplitter(sfInfo);

    //theSplitter.split(samples, correlator.getCorrelationZeroes(), instrument);

    //Splitter theSplitter(correlator, startSample);

	return 0;
}