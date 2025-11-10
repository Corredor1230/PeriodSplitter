#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include <iostream>
#include <string>
#include <thread>
#include <stdexcept>

#include "gui/EditorViewModel.h"
#include "gui/SihatEditor.h"
#include "file/SihatJson.h" // For loading/saving

// --- This is your DSP logic, running on its own thread ---
// (We keep it separate, often in its own file)
void runAnalysisThread(
    bool isBulk,
    std::string path,
    Sitrano::AnalysisConfig config, // Passed by-value (COPY)
    Sitrano::Settings settings,     // Passed by-value (COPY)
    SihatFile::OutInfo info,        // Passed by-value (COPY)
    std::atomic<bool>& isProcessingFlag)
{
    std::cout << "Analysis thread started..." << std::endl;
    try
    {
        if (isBulk)
        {
            SihatFile::processFolder(path, info, config, settings);
        }
        else
        {
            SihatFile::processFile(path, info, config, settings);
        }
        std::cout << "Analysis thread finished." << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Analysis thread failed: " << e.what() << std::endl;
    }

    // Signal that processing is done
    isProcessingFlag = false;
}


// --- Main Application ---
int main()
{
    // --- 1. Setup ---
    EditorViewModel viewModel; // The one and only "bridge" object
    viewModel.jsonPath = "../source/configFiles/default.json"; // Or get from argv

    // Load initial settings into the view model
    try
    {
        SihatJson::loadSettings(viewModel.config, viewModel.settings, viewModel.info, viewModel.jsonPath);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Failed to load settings: " << e.what() << ". Using defaults." << std::endl;
        // ... you might set defaults here ...
    }

    // Initialize GLFW & ImGui
    GLFWwindow* window = SihatEditor::initializeGlfw("Sitrano Editor", 1280, 720);
    if (!window) return -1;
    SihatEditor::initializeImgui(window);

    std::thread analysisThread; // The processing thread
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // --- 2. Wire Up Callbacks (This is the Controller Logic) ---
    // We "teach" the view model what to do when buttons are pressed.

    viewModel.onSaveRequested = [&]() {
        std::cout << "Save requested..." << std::endl;
        SihatJson::saveSettings(viewModel.config, viewModel.settings, viewModel.info, viewModel.jsonPath);
        std::cout << "Settings saved to " << viewModel.jsonPath << std::endl;
        };

    // Helper lambda to launch analysis
    auto run_logic = [&](bool isBulk) {
        if (viewModel.isProcessing.load()) return; // Don't start a new job

        // File dialogs MUST run on the main thread
        std::string path = isBulk ? SihatFile::openFolderDialog() : SihatFile::openFileDialog();
        if (path.empty())
        {
            std::cout << "Analysis cancelled (no file/folder selected)." << std::endl;
            return;
        }

        std::cout << "Starting analysis on: " << path << std::endl;
        viewModel.isProcessing = true; // Set flag so GUI updates

        // Clean up old thread if it's finished
        if (analysisThread.joinable())
            analysisThread.join();

        // Launch new thread with COPIES of data
        analysisThread = std::thread(
            runAnalysisThread,
            isBulk,
            path,
            viewModel.config,     // Copied
            viewModel.settings,   // Copied
            viewModel.info,       // Copied
            std::ref(viewModel.isProcessing) // Pass flag by reference
        );
        };

    viewModel.onRunSingleRequested = [&]() { run_logic(false); };
    viewModel.onRunBulkRequested = [&]() { run_logic(true); };


    // --- 3. Main Application Loop ---
    while (!glfwWindowShouldClose(window))
    {
        // Poll for and process events (mouse, keyboard)
        glfwPollEvents();

        // Start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- RENDER THE VIEW ---
        // This one function call draws the entire GUI.
        SihatEditor::Render(viewModel);

        // --- Rendering ---
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);

        // --- Post-render logic ---
        // Check if the thread finished
        if (!viewModel.isProcessing && analysisThread.joinable())
        {
            analysisThread.join();
            std::cout << "Analysis thread joined." << std::endl;
        }
    }

    // --- 4. Cleanup ---
    if (analysisThread.joinable())
    {
        analysisThread.join(); // Wait for thread to finish
    }
    SihatEditor::cleanup(window);
    return 0;
}