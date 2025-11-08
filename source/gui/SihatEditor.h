#pragma once

#include<iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

// --- ImGui ---
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "include/SitranoHeader.h"
#include "file/SihatFile.h"

class SihatEditor
{
public:

	SihatEditor() {};
	~SihatEditor() {};

	static void glfwErrorCallback(int error, const char* description)
	{
		std::cerr << "GLFW error" << ": " << description << std::endl;
	}

	GLFWwindow* initializeGlfw(const char* title, int width, int height);

	void initializeImgui(GLFWwindow* window);

	void renderLoop(GLFWwindow* window, Sitrano::AnalysisConfig& config,
		Sitrano::Settings& settings, SihatFile::OutInfo& info,
		const std::string& jsonPath);

	void cleanup(GLFWwindow* window);
};