#include "SihatEditor.h"

#include <string>
#include <stdexcept>
#include "file/SihatJson.h"

#include "ConfigEditor.h"

GLFWwindow* SihatEditor::initializeGlfw(const char* title, int width, int height)
{
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return nullptr;
    }

    // Set GL version (ImGui works well with 3.3)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // --- Initialize GLAD ---
    // This MUST be done after making the context current
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        glfwDestroyWindow(window);
        glfwTerminate();
        return nullptr;
    }

    return window;
}

void SihatEditor::initializeImgui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;    
    //io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;   

    // Setup ImGui style
    ImGui::StyleColorsDark();

    // When viewports are enabled, tweak WindowRounding/WindowBg so platform windows can look identical
    /*ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }*/

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
}

void SihatEditor::renderLoop(GLFWwindow* window,
    Sitrano::AnalysisConfig& config,
    Sitrano::Settings& settings,
    SihatFile::OutInfo& info,
    const std::string& jsonPath)
{
    bool bulkProcess = config.bulkProcess;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    while (!glfwWindowShouldClose(window))
    {
        // Poll for and process events (mouse, keyboard)
        glfwPollEvents();

        // Start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- Main GUI Window ---
        {
            // We set a min size for the window
            ImGui::SetNextWindowSize(ImVec2(600, 400), ImGuiCond_FirstUseEver);
            ImGui::Begin("Sitrano Analysis Configuration");

            // Display the path to the loaded JSON
            ImGui::Text("Settings loaded from: %s", jsonPath.c_str());
            ImGui::Separator();

            // --- Control Panel Buttons ---
            if (ImGui::Button("Save Settings"))
            {
                SihatJson::saveSettings(config, settings, info, jsonPath);
                std::cout << "Settings saved to " << jsonPath << std::endl;
            }
            ImGui::SameLine();
            ImGui::Checkbox("Bulk Process Folder", &bulkProcess);
            ImGui::SameLine();

            // --- THIS IS THE LOGIC FROM THE OLD MAIN ---
            if (ImGui::Button("Run Analysis"))
            {
                std::cout << "Run Analysis Clicked..." << std::endl;
                if (!bulkProcess)
                {
                    std::string filename = SihatFile::openFileDialog();
                    if (!filename.empty())
                    {
                        std::cout << "Processing single file: " << filename << std::endl;
                        SihatFile::processFile(filename, info, config, settings);
                    }
                }
                else
                {
                    std::string inputDir = SihatFile::openFolderDialog();
                    // Your old main had this, you can uncomment if needed
                    if (!inputDir.empty())
                    {
                        std::cout << "Processing folder: " << inputDir << std::endl;
                        SihatFile::processFolder(inputDir, info, config, settings);
                    }
                }
            }
            // --- END OF OLD MAIN LOGIC ---

            ImGui::Separator();

            // --- Settings Editor ---
            // Create a child window to make the settings scrollable
            ImGui::BeginChild("SettingsChildWindow", ImVec2(0, 0), true);

            // Call our editor class to draw the contents!
            ConfigEditor::Render(config, settings, info);

            ImGui::EndChild();

            ImGui::End();
        }

        // --- Rendering ---
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }
}

void SihatEditor::cleanup(GLFWwindow* window)
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}