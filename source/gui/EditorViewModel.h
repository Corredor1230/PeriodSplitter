#pragma once

#include "include/SitranoHeader.h" 
#include "include/ResynthConfig.h"
#include "include/InterpConfig.h"
#include "file/SihatFile.h"   
#include "support/SihatLogger.h"
#include <string>
#include <atomic>     // For the isProcessing flag
#include <functional> // For the callbacks

enum class AppMode
{
    Analysis,
    Resynthesis,
    Interpolation
};

struct EditorViewModel
{

    AppMode currentMode = AppMode::Analysis;

    // The GUI will read/write these values directly.
    Sihat::AnalysisConfig   config;
    ResynthConfig           rConfig;
    InterpConfig::Config iConfig;
    InterpConfig::Parameters iParams;
    Sihat::Settings         settings;
    SihatFile::OutInfo      info;
    std::string             jsonPath;

    // The GUI reads this to disable buttons, show spinners, etc.
    std::atomic<bool> isProcessing{ false };

    bool hasAudioToPlay = false;
    std::vector<float> generatedAudio;

    // The GUI will call these functions when buttons are clicked.
    // Our 'SihatApplication.cpp' will provide the definitions for these.
    std::function<void(void)> onSaveRequested;
    std::function<void(void)> onRunSingleRequested;
    std::function<void(void)> onRunBulkRequested;
    std::function<void(void)> onRunResynthRequested;

    std::function<void(void)> onLoadFileA;
    std::function<void(void)> onLoadFileB;
    std::function<void(void)> onAlphaChanged;
    std::function<void(void)> onPlayAudio;
    std::function<void(void)> onSaveAudio;

    SihatLogger logger;
};