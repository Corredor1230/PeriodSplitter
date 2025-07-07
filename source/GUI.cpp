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
    if (glfwGetKey(m_window.get(), GLFW_KEY_RIGHT) == GLFW_PRESS && m_currentWindow < m_windows.size() - 1)
        m_currentWindow++;
    if (glfwGetKey(m_window.get(), GLFW_KEY_LEFT) == GLFW_PRESS && m_currentWindow > 0)
        m_currentWindow--;
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
