#pragma once

#include "glad/glad.h" // Or your GL loader
#include "GLFW/glfw3.h"
#include "EditorViewModel.h" // Include the new view model

class SihatEditor
{
public:
    /**
     * @brief Initializes GLFW and GLAD, creates a window.
     */
    static GLFWwindow* initializeGlfw(const char* title, int width, int height);

    /**
     * @brief Initializes ImGui for an OpenGL/GLFW context.
     */
    static void initializeImgui(GLFWwindow* window);

    /**
     * @brief Renders the editor GUI for one frame.
     * Call this inside your main loop, between ImGui::NewFrame() and ImGui::Render().
     *
     * @param viewModel The application's view model.
     */
    static void Render(EditorViewModel& viewModel);

    /**
     * @brief Shuts down ImGui and GLFW.
     */
    static void cleanup(GLFWwindow* window);
};