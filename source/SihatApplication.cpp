#include "SihatApplication.h"

// All the includes that main.cpp used to have
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "include/ResynthHeader.h"
#include "include/ResynthConfig.h"
#include "resynthesis/SihatResynth.h"
#include "include/InterpConfig.h"
#include "include/InterpHeader.h"
#include "interpolation/Interpolator.h"

#include <iostream>
#include <stdexcept>
#include "file/SihatJson.h"
#include "gui/SihatEditor.h"
#include "support/SihatLogger.h"

void runAnalysisThread(
    bool isBulk,
    std::string path,
    Sihat::AnalysisConfig config,
    Sihat::Settings settings,
    SihatFile::OutInfo info,
    std::atomic<bool>& isProcessingFlag,
    SihatLogger& logger)
{
    logger.logTemp("Analysis started...");
    logger.log("Processing: " + path);
    try
    {
        if (isBulk) {
            SihatFile::processFolder(path, info, config, settings);
            logger.log("Folder processing complete.");
        }
        else {
            SihatFile::processFile(path, info, config, settings);
            logger.log("File processing complete.");
        }
        logger.logTemp("Analysis finished successfully.");
    }
    catch(const std::exception& e)
    {
        logger.log("[ERROR] " + std::string(e.what()));
        logger.logTemp("Analysis failed. Check log for details.");
    }
    isProcessingFlag = false;
}

void runResynthThread(
    std::string path,
    std::atomic<bool>& isProcessingFlag,
    SihatLogger& logger,
    const ResynthConfig& config)
{
    logger.logTemp("Resynthesis started...");
    logger.log("Processing: " + path);
    
    Synth::Sihat sihat;   

    try
    {
        sihat = SihatFile::loadSihatFile(path);
        Resynthesizer synth(sihat, config);
        synth.resynthesize();
        logger.log("Resynthesis complete.");
    }
    catch(const std::exception& e)
    {
        logger.log("[ERROR] " + std::string(e.what()));
        logger.logTemp("Resynthesis failed. Check log for details.");
    }

    isProcessingFlag = false;
    
}

void runInterpolationThread(
    std::string pathA,
    std::string pathB,
    InterpConfig::Config& iConfig,
    InterpConfig::Parameters& iParams,
    float alpha,
    std::atomic<bool>& isProcessingFlag,
    SihatLogger& logger,
    bool& hasAudioFlag,
    std::vector<float>& outAudio) // Passed by ref so the GUI can access it later
{
    logger.logTemp("Interpolating and Resynthesizing...");
    
    try {
        // 1. Load the files
        Synth::Sihat sihatA = SihatFile::loadSihatFile(pathA);
        Synth::Sihat sihatB = SihatFile::loadSihatFile(pathB);

        // 2. Interpolate (using the aggregator we built earlier)
        Interpolator interp(iConfig, iParams);
        // interp.conf.alpha = alpha; // Set up your conf however you structured it
        float targetF0 = iConfig.f0; // Example
        Synth::Sihat interpolatedSihat = interp.interpolateSihat(sihatA, sihatB, targetF0);

        // 3. Resynthesize
        Resynthesizer synth(interpolatedSihat, ResynthConfig()); // Pass relevant config
        synth.resynthesize(); // Assuming this generates the audio internally
        
        // 4. Extract the audio
        // You'll need a method in your Resynthesizer to get the resulting buffer
        outAudio = synth.getResynthAudio();
        
        hasAudioFlag = true;
        logger.log("Interpolation and audio generation complete.");
    }
    catch(const std::exception& e) {
        logger.log("[ERROR] Interpolation failed: " + std::string(e.what()));
        hasAudioFlag = false;
    }

    isProcessingFlag = false;
}

SihatApplication::SihatApplication(const std::string& jsonPath)
{
    viewModel.jsonPath = jsonPath;

    try {
        SihatJson::loadSettings(viewModel.config, viewModel.settings,
            viewModel.info, viewModel.rConfig, viewModel.jsonPath);
    }
    catch (const std::exception& e) {
        //std::cerr << "Failed to load settings: " << e.what() << ". Using defualts." << std::endl;
        viewModel.logger.log("[ERROR] Failed to load settings: " + std::string(e.what()));
    }

    m_audioPlayer.init(viewModel.iConfig.sampleRate);

    initialize();

    viewModel.onSaveRequested = [this]() {
        this->onSave();
        };

    viewModel.onRunSingleRequested = [this]() {
        this->onRun(false);
        };

    viewModel.onRunBulkRequested = [this]() {
        this->onRun(true);
        };

    viewModel.onRunResynthRequested = [this]() {
        this->onRunResynth();
        };
    viewModel.onLoadFileA = [this]() {
        std::string path = SihatFile::openSihatDialog();
        if (!path.empty()) viewModel.iConfig.pathA = path;
        };

    viewModel.onLoadFileB = [this]() {
        std::string path = SihatFile::openSihatDialog();
        if (!path.empty()) viewModel.iConfig.pathB = path;
        };

    viewModel.onAlphaChanged = [this]() {
        this->onRunInterpolation();
        };

    viewModel.onPlayAudio = [this]() {
        viewModel.logger.logTemp("Playing audio...");
        m_audioPlayer.playAudio(viewModel.generatedAudio);
        };

    viewModel.onSaveAudio = [this]() {
        std::string savePath = SihatFile::openFolderDialog(); // Implement this if you haven't
        if (!savePath.empty()) {
            // Example: AudioFile::saveWav(savePath, viewModel.generatedAudio, sampleRate);
            viewModel.logger.log("Audio saved to " + savePath);
        }
        };
}

SihatApplication::~SihatApplication() {
    shutdown();
}

void SihatApplication::initialize()
{
    m_window = SihatEditor::initializeGlfw("sihat Editor", 1280, 720);
    if (!m_window) {
        throw std::runtime_error("Failed to initialize Glfw");
    }
    SihatEditor::initializeImgui(m_window);
    viewModel.logger.logTemp("Ready.");
}

void SihatApplication::shutdown()
{
    checkAndJoinThread();
    SihatEditor::cleanup(m_window);
}

void SihatApplication::run()
{
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    while (!glfwWindowShouldClose(m_window))
    {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        SihatEditor::Render(viewModel);
        ImGui::Render();

        int displayW, displayH;
        glfwGetFramebufferSize(m_window, &displayW, &displayH);
        glViewport(0, 0, displayW, displayH);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(m_window);

        checkAndJoinThread();
    }
}

void SihatApplication::checkAndJoinThread()
{
    if (!viewModel.isProcessing && analysisThread.joinable())
    {
        analysisThread.join();
        viewModel.logger.log("Analysis thread joined");
    }
}

void SihatApplication::onSave()
{
    viewModel.logger.logTemp("Saving settings...");
    try {
        SihatJson::saveSettings(viewModel.config, viewModel.settings, viewModel.info, viewModel.rConfig, viewModel.jsonPath);
        viewModel.logger.log("Settings saved to " + viewModel.jsonPath);
        viewModel.logger.logTemp("Ready.");
    }
    catch (const std::exception& e)
    {
        viewModel.logger.log("[ERROR] Failed to save settings: " + std::string(e.what()));
        viewModel.logger.logTemp("Save failed.");
    }
}

void SihatApplication::onRun(bool isBulk)
{
    if (viewModel.isProcessing.load()) return;

    std::string path = isBulk ? SihatFile::openFolderDialog() : SihatFile::openFileDialog();
    if (path.empty()) {
        viewModel.logger.logTemp("Analysis cancelled");
        return;
    }

    viewModel.logger.logTemp("Starting analysis on " + path);
    viewModel.isProcessing = true;

    checkAndJoinThread();

    analysisThread = std::thread(
        runAnalysisThread,
        isBulk,
        path,
        viewModel.config,
        viewModel.settings,
        viewModel.info,
        std::ref(viewModel.isProcessing),
        std::ref(viewModel.logger)
    );
}

void SihatApplication::onRunResynth()
{
    if (viewModel.isProcessing.load()) return;

    std::string path = SihatFile::openSihatDialog();
    if (path.empty()) {
        viewModel.logger.logTemp("Resynthesis cancelled");
        return;
    }

    viewModel.logger.logTemp("Starting resynthesis on " + path);
    viewModel.isProcessing = true;

    checkAndJoinThread();

    analysisThread = std::thread(
        runResynthThread,
        path,
        std::ref(viewModel.isProcessing),
        std::ref(viewModel.logger),
        viewModel.rConfig
    );
}

void SihatApplication::onRunInterpolation()
{
    if (viewModel.isProcessing.load()) return;
    
    // Reset audio state
    viewModel.hasAudioToPlay = false;
    viewModel.generatedAudio.clear();
    viewModel.isProcessing = true;
    
    checkAndJoinThread();

    analysisThread = std::thread(
        runInterpolationThread,
        viewModel.iConfig.pathA,
        viewModel.iConfig.pathB,
        std::ref(viewModel.iConfig),
        std::ref(viewModel.iParams),
        viewModel.iConfig.alpha,
        std::ref(viewModel.isProcessing),
        std::ref(viewModel.logger),
        std::ref(viewModel.hasAudioToPlay),
        std::ref(viewModel.generatedAudio)
    );
}