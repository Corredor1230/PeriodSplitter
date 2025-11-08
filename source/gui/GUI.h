#pragma once

#include <vector>
#include <utility>
#include <string>
#include <memory>
#include <GLFW/glfw3.h>

class WaveformViewer {
public:
    WaveformViewer(const std::vector<float>& samples, const std::vector<std::pair<int, int>>& windows);
    ~WaveformViewer();
    void run(); // Starts GUI loop

private:

    double m_keyPressTimeLeft = 0.0;
    double m_keyPressTimeRight = 0.0;
    double m_nextChangeTimeLeft = 0.0;
    double m_nextChangeTimeRight = 0.0;

    void initWindow();
    void initImGui();
    void handleInput();
    void renderWaveform();
    void cleanup();

    std::vector<float> m_samples;
    std::vector<std::pair<int, int>> m_windows;
    int m_currentWindow = 0;

    using GLFWwindowPtr = std::unique_ptr<GLFWwindow, decltype(&glfwDestroyWindow)>;
    GLFWwindowPtr m_window;
};
