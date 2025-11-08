#include<iostream>
#include<string.h>
#include<vector>
#include"file/SihatFile.h"
#include"file/JsonAdapters.h"
#include"file/SihatJson.h"

int main()
{
    bool bulkProcess = true;

    Sitrano::AnalysisConfig config;
    Sitrano::Settings settings;
    SihatFile::OutInfo info;

    std::string jsonPath = SihatFile::openJsonDialog();
    SihatJson::loadSettings(config, settings, info, jsonPath);
    //SihatJson::saveSettings(config, settings, info, jsonPath);


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
