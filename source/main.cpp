#include<iostream>
#include<string.h>
#include<vector>
#include"file/SihatFile.h"
#include"file/JsonAdapters.h"
#include"file/SihatJson.h"
#include"gui/SihatEditor.h"

int main()
{
    glfwSetErrorCallback(SihatEditor::glfwErrorCallback);

    // --- Load Data FIRST ---
    // This logic is from your old main()
    bool useDefault = true; // You can change this
    Sitrano::AnalysisConfig config;
    Sitrano::Settings settings;
    SihatFile::OutInfo info;
    std::string jsonPath;

    if (useDefault)
    {
        jsonPath = "C:/Users/usuario/Documents/Programming/CMake_Learning/PeriodicSplitter/source/configFiles/default.json";
    }
    else
    {
        jsonPath = SihatFile::openJsonDialog();
    }

    if (jsonPath.empty())
    {
        std::cerr << "No JSON file selected. Exiting." << std::endl;
        return -1;
    }
    SihatJson::loadSettings(config, settings, info, jsonPath);
    std::cout << "Loaded settings from " << jsonPath << std::endl;

    // --- Create an instance of the Editor ---
    SihatEditor editor;

    // --- Initialize Window & GUI ---
    GLFWwindow* window = editor.initializeGlfw("Sitrano Analysis", 1280, 720);
    if (!window)
    {
        return -1; // Error already printed
    }

    editor.initializeImgui(window);

    // --- Run Application ---
    // This one function now runs the whole app.
    // It will loop until the window is closed.
    editor.renderLoop(window, config, settings, info, jsonPath);

    // --- Cleanup ---
    editor.cleanup(window);

    return 0;
}
