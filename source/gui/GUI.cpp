#include "GUI.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdexcept>
#include <iostream>

WaveformViewer::WaveformViewer(const std::vector<float>& samples,
    const std::vector<std::pair<int, int>>& windows)
    : m_samples(samples), m_windows(windows),
    m_window(nullptr, glfwDestroyWindow)
{
    initWindow();
    initImGui();
}

WaveformViewer::~WaveformViewer() {
    cleanup();
}

void WaveformViewer::initWindow() {
    if (!glfwInit())
        throw std::runtime_error("Failed to initialize GLFW");

    GLFWwindow* raw = glfwCreateWindow(800, 600, "Waveform Viewer", nullptr, nullptr);
    if (!raw)
        throw std::runtime_error("Failed to create GLFW window");

    m_window.reset(raw);
    glfwMakeContextCurrent(m_window.get());
    glfwSwapInterval(1); // Enable vsync
}

void WaveformViewer::initImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(m_window.get(), true);
    ImGui_ImplOpenGL3_Init("#version 130");
}

void WaveformViewer::handleInput() {
    // --- Configuration for key repeat behavior ---
    const double INITIAL_DELAY = 0.4;       // seconds to wait before repeating starts
    const double ACCELERATION_TIME = 1.0;     // seconds to hold before speed increases
    const double REPEAT_RATE_SLOW = 0.1;      // seconds between changes (the slow speed)
    const double REPEAT_RATE_FAST = 0.03;     // seconds between changes (the fast speed)

    double currentTime = glfwGetTime();

    // --- Handle Right Arrow Key ---
    // A. This block handles the very first press
    if (ImGui::IsKeyPressed(ImGuiKey_RightArrow, false)) {
        if (m_currentWindow < m_windows.size() - 1) {
            m_currentWindow++;
        }
        // Record the press time and schedule the next change after the initial delay
        m_keyPressTimeRight = currentTime;
        m_nextChangeTimeRight = currentTime + INITIAL_DELAY;
    }
    // B. This block handles a sustained press (holding the key down)
    else if (ImGui::IsKeyDown(ImGuiKey_RightArrow)) {
        if (currentTime >= m_nextChangeTimeRight) {
            if (m_currentWindow < m_windows.size() - 1) {
                m_currentWindow++;
            }
            // Check how long the key has been held to determine the repeat rate
            double holdDuration = currentTime - m_keyPressTimeRight;
            double currentRate = (holdDuration > ACCELERATION_TIME) ? REPEAT_RATE_FAST : REPEAT_RATE_SLOW;

            // Schedule the next change
            m_nextChangeTimeRight = currentTime + currentRate;
        }
    }

    // --- Handle Left Arrow Key (same logic) ---
    if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow, false)) {
        if (m_currentWindow > 0) {
            m_currentWindow--;
        }
        m_keyPressTimeLeft = currentTime;
        m_nextChangeTimeLeft = currentTime + INITIAL_DELAY;
    }
    else if (ImGui::IsKeyDown(ImGuiKey_LeftArrow)) {
        if (currentTime >= m_nextChangeTimeLeft) {
            if (m_currentWindow > 0) {
                m_currentWindow--;
            }
            double holdDuration = currentTime - m_keyPressTimeLeft;
            double currentRate = (holdDuration > ACCELERATION_TIME) ? REPEAT_RATE_FAST : REPEAT_RATE_SLOW;
            m_nextChangeTimeLeft = currentTime + currentRate;
        }
    }
}

void WaveformViewer::renderWaveform() {
    ImGui::Begin("Waveform Viewer");

    auto [start, end] = m_windows[m_currentWindow];
    ImGui::Text("Window %d: [%d - %d]", m_currentWindow, start, end);

    if (end > start && end < m_samples.size()) {
        ImGui::PlotLines("Waveform", m_samples.data() + start, end - start, 0, nullptr, -1.0f, 1.0f, ImVec2(700, 200));
    }

    ImGui::End();
}

void WaveformViewer::cleanup() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(m_window.release());
    glfwTerminate();
}

void WaveformViewer::run() {
    while (!glfwWindowShouldClose(m_window.get())) {
        glfwPollEvents();
        handleInput();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        renderWaveform();

        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(m_window.get(), &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(m_window.get());
    }
}
