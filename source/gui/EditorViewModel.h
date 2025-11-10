#pragma once

#include "include/SitranoHeader.h" 
#include "file/SihatFile.h"   
#include "support/SihatLogger.h"
#include <string>
#include <atomic>     // For the isProcessing flag
#include <functional> // For the callbacks

struct EditorViewModel
{
    // The GUI will read/write these values directly.
    Sitrano::AnalysisConfig config;
    Sitrano::Settings       settings;
    SihatFile::OutInfo      info;
    std::string             jsonPath;

    // The GUI reads this to disable buttons, show spinners, etc.
    std::atomic<bool> isProcessing{ false };

    // The GUI will call these functions when buttons are clicked.
    // Our 'main.cpp' will provide the definitions for these.
    std::function<void(void)> onSaveRequested;
    std::function<void(void)> onRunSingleRequested;
    std::function<void(void)> onRunBulkRequested;

    SihatLogger logger;
};