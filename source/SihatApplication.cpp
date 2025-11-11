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

SihatApplication::SihatApplication(const std::string& jsonPath)
{
    viewModel.jsonPath = jsonPath;

    try {
        SihatJson::loadSettings(viewModel.config, viewModel.settings,
            viewModel.info, viewModel.jsonPath);
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load settings: " << e.what() << ". Using defualts." << std::endl;
        viewModel.logger.log("[ERROR] Failed to load settings: " + std::string(e.what()));
    }

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
        SihatJson::saveSettings(viewModel.config, viewModel.settings, viewModel.info, viewModel.jsonPath);
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