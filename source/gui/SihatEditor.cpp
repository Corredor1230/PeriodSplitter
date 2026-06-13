#include "SihatEditor.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"  
#include "imgui_impl_opengl3.h" 
#include "ConfigEditor.h"
#include "imgui_internal.h"

static void renderLogger(SihatLogger& logger)
{
    {
        std::lock_guard<std::mutex> lock(logger.getMutex());
        ImGui::Text("Status: %s", logger.getTemporaryMessage().c_str());
    }

    ImGui::Separator();

    if (ImGui::CollapsingHeader("Analysis log"))
    {
        float logHeight = 200.f;
        ImGui::BeginChild("Permanent log", ImVec2(0, logHeight), true, ImGuiWindowFlags_HorizontalScrollbar);

        std::lock_guard<std::mutex> lock(logger.getMutex());
        const auto& messages = logger.getPermanentMessages();

        ImGuiListClipper clipper;
        clipper.Begin(messages.size());
        while (clipper.Step())
        {
            for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; i++)
            {
                ImGui::TextUnformatted(messages[i].c_str());
            }
        }
        clipper.End();

        if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
        {
            ImGui::SetScrollHereY(1.0f);
        }

        ImGui::EndChild();
    }
}

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
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  

    // Setup ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
}

void SihatEditor::Render(EditorViewModel& viewModel)
{
    // --- Main GUI Window ---
    ImGui::SetNextWindowSize(ImVec2(600, 400), ImGuiCond_FirstUseEver);
    ImGui::Begin("Sihat Configuration");
    {
        ImGui::Text("Settings loaded from: %s", viewModel.jsonPath.c_str());
        ImGui::Separator();

        // --- Control Panel ---
        bool processing = viewModel.isProcessing.load();

        // Disable buttons if processing
        ImGui::BeginDisabled(processing);
        if (ImGui::BeginTabBar("AppModeTabs"))
        {
            if (ImGui::BeginTabItem("Analysis"))
            {
                viewModel.currentMode = AppMode::Analysis;
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Resynthesis"))
            {
                viewModel.currentMode = AppMode::Resynthesis;
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Interpolation"))
            {
                viewModel.currentMode = AppMode::Interpolation;
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::EndDisabled();
        ImGui::Separator();

        ImGui::BeginDisabled(processing);

        if (ImGui::Button("Save Settings"))
        {
            // 1. Just call the callback
            if (viewModel.onSaveRequested)
                viewModel.onSaveRequested();
        }

        ImGui::SameLine();
        
        if (viewModel.currentMode == AppMode::Analysis)
        {
            ImGui::Checkbox("Bulk Process Folder", &viewModel.config.bulkProcess);

            ImGui::SameLine();
            if (ImGui::Button(processing ? "Processing..." : "Run Analysis"))
            {
                // 3. Call the correct callback based on state
                if (viewModel.config.bulkProcess)
                {
                    if (viewModel.onRunBulkRequested)
                        viewModel.onRunBulkRequested();
                }
                else
                {
                    if (viewModel.onRunSingleRequested)
                        viewModel.onRunSingleRequested();
                }
            }
        }
        else if (viewModel.currentMode == AppMode::Resynthesis)
        {
            if (ImGui::Button(processing ? "Processing Resynth..." : "Run Resynthesis"))
            {
                if (viewModel.onRunResynthRequested) viewModel.onRunResynthRequested();
            }
        }

        ImGui::EndDisabled(); // Re-enable GUI

        ImGui::Separator();

        if (processing)
        {
            ImGui::Text(viewModel.currentMode == AppMode::Analysis ? "Analysis in progress, please wait..." : "Resynthesis in progress, please wait...");
            // You could add an ImGui::Spinner or loading bar here
        }

        // --- Settings Editor ---
        ImGui::BeginChild("SettingsChildWindow", ImVec2(0, 0), true);
        {
            renderLogger(viewModel.logger);
            ImGui::Separator();

            if (viewModel.currentMode == AppMode::Analysis)
            {
                ConfigEditor::Render(viewModel.config, viewModel.settings, viewModel.info);
            }
            else if (viewModel.currentMode == AppMode::Resynthesis)
            {
                ConfigEditor::ResynthRender(viewModel.rConfig);
            }
            else if (viewModel.currentMode == AppMode::Interpolation)
            {
                ConfigEditor::InterpolationRender(viewModel);
            }
        }

        ImGui::EndChild();
    }
    ImGui::End();
}

void SihatEditor::cleanup(GLFWwindow* window)
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}