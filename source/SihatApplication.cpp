#include "SihatApplication.h"

// All the includes that main.cpp used to have
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include <iostream>
#include <stdexcept>
#include "file/SihatJson.h"
#include "gui/SihatEditor.h"
#include "support/SihatLogger.h"

void runAnalysisThread(
    bool isBulk,
    std::string path,
    Sitrano::AnalysisConfig config,
    Sitrano::Settings settings,
    SihatFile::OutInfo info,
    std::atomic<bool>& isProcessingFlag,
    SihatLogger& logger)
{

}