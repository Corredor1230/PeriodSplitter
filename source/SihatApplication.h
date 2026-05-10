#pragma once

#include"gui/EditorViewModel.h"
#include<thread>
#include<string>

struct GLFWwindow;

class SihatApplication
{
public:
	SihatApplication(const std::string& jsonPath);
	~SihatApplication();
	void run();

private:
	void initialize();
	void shutdown();
	void checkAndJoinThread();
	void onSave();
	void onRun(bool isBulk);
	void onRunResynth();

	GLFWwindow* m_window = nullptr;
	EditorViewModel viewModel;
	std::thread analysisThread;
};