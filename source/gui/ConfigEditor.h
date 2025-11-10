#pragma once 

#include"imgui.h"
#include"imgui_stdlib.h"
#include<string>
#include"include/SitranoHeader.h"
#include"file/SihatFile.h"
#include"file/SihatJson.h"


class ConfigEditor
{
public:
    /**
     * @brief Renders the configuration editor window.
     * Call this inside your main GLFW loop, between ImGui::NewFrame() and ImGui::Render().
     *
     * @param config Reference to your main config struct.
     * @param settings Reference to your settings struct.
     * @param info Reference to your file info struct.
     * @param p_open Pointer to a bool to control the window's visibility (optional).
     */
    static void Render(Sitrano::AnalysisConfig& config, Sitrano::Settings& settings, SihatFile::OutInfo& info, bool* p_open = nullptr)
    {
        // Begin the ImGui window
        if (!ImGui::Begin("Analysis Configuration Editor", p_open))
        {
            // Window is collapsed, so we stop drawing anything inside it.
            ImGui::End();
            return;
        }

        // --- 1. Analysis Modules (Settings) ---
        if (ImGui::CollapsingHeader("Analysis Modules", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Pitch Analysis", &settings.pitchAnalysis);
            ImGui::Checkbox("Transient Separation", &settings.transientSeparation);
            ImGui::Checkbox("Period Analysis", &settings.periodAnalysis);
            ImGui::Checkbox("Overtone Analysis", &settings.overtoneAnalysis);
            ImGui::Checkbox("Harmonic Analysis", &settings.harmonicAnalysis);
            ImGui::Checkbox("Noise Analysis", &settings.noiseAnalysis);
        }

        // --- 2. Output File Info (Info) ---
        if (ImGui::CollapsingHeader("Output File Info", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Uses imgui_stdlib.h to allow ImGui::InputText to work with std::string
            ImGui::InputText("Output Directory", &info.outDir);
            ImGui::InputText("File Prefix", &info.prefix);
            ImGui::InputText("File Extension", &info.extension);
        }

        // --- 3. Main Configuration (Config) ---
        if (ImGui::CollapsingHeader("Main Configuration", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // --- Basic Parameters ---
            ImGui::InputInt("N (FFT Size)", &config.nfft);
            ImGui::InputInt("Hop Size", &config.hopSize);
            ImGui::InputInt("Start Sample", &config.startSample);
            ImGui::InputInt("Max Harmonics", &config.numHarmonics);
            ImGui::InputFloat("Global Tolerance", &config.tolerance);
            ImGui::Checkbox("Verbose Logging", &config.verbose);
            ImGui::Checkbox("Bulk Process", &config.bulkProcess);

            ImGui::Separator();

            // --- Nested Settings Groups ---

            if (ImGui::TreeNode("Correlation Settings"))
            {
                ImGui::InputFloat("Offset", &config.cSettings.periodStartOffsetMs);
                ImGui::InputFloat("Threshold", &config.cSettings.correlationThreshold);
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Harmonic Settings"))
            {
                ImGui::Checkbox("Enabled", &config.hSettings.applyHanning);
                ImGui::InputFloat("Tolerance", &config.hSettings.toleranceValue);
                const char* windowStyles[] = { "Period Loop", "Single Period", "Audio Chunk" };
                int styleCount = IM_ARRAYSIZE(windowStyles);

                int styleAsInt = static_cast<int>(config.hSettings.style);
                if (ImGui::Combo("Window Style", &styleAsInt, windowStyles, styleCount))
                {
                    config.hSettings.style = static_cast<Sitrano::WindowStyle>(styleAsInt);
                }

                // I'm adding a small tooltip to remind you to change the names
                ImGui::SameLine();
                ImGui::TextDisabled("(?)");
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
                {
                    ImGui::BeginTooltip();
                    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
                    ImGui::TextUnformatted("NOTE: I've assumed the names and order for the 'Window Style' enum (e.g., 0=Rectangular, 1=Hann, ...).\n\nYou'll need to update the 'windowStyles' array in ConfigEditor.h to match your actual enum definition!");
                    ImGui::PopTextWrapPos();
                    ImGui::EndTooltip();
                }
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Overtone Settings"))
            {
                ImGui::InputInt("First Sample", &config.oSettings.overtoneFirstSample);
                ImGui::InputInt("MFFT", &config.oSettings.fftSize);
                ImGui::InputFloat("Ignore Threshold", &config.oSettings.overtoneThreshold);
                ImGui::InputDouble("Tolerance", &config.oSettings.toleranceValue);
                ImGui::Checkbox("Use Tolerance", &config.oSettings.useTolerance);
                ImGui::Checkbox("Absolute Threshold", &config.oSettings.setAbsoluteThreshold);
                ImGui::Checkbox("Post Transient Start", &config.oSettings.postTransientStart);
                ImGui::Checkbox("Sum Amplitudes", &config.oSettings.sumAmplitudes);
                ImGui::Checkbox("Use Custom Signal", &config.oSettings.useCustomSignal);
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Pitch Settings"))
            {
                ImGui::InputFloat("Min Frequency", &config.pSettings.minFreq);
                ImGui::InputFloat("Max Frequency", &config.pSettings.maxFreq);
                ImGui::InputFloat("Threshold", &config.pSettings.modeThreshold);
                ImGui::InputFloat("In Cents", &config.pSettings.toleranceInCents);
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Transient Settings"))
            {
                ImGui::InputInt("Start Sample", &config.tSettings.tStartSample);
                ImGui::InputFloat("Threshold", &config.tSettings.transientThreshold);
                ImGui::InputFloat("Factor", &config.tSettings.transientFactor);
                ImGui::InputInt("Pre-Attack", &config.tSettings.preAttack);
                ImGui::InputInt("RMS Sample Size", &config.tSettings.rmsSampleSize);
                ImGui::InputInt("RMS Hop Length", &config.tSettings.rmsSampleHopLength);
                ImGui::InputFloat("RMS Hop Ratio", &config.tSettings.transientRmsHopRatio);
                ImGui::InputFloat("RMS Size (ms)", &config.tSettings.transientRmsSizeMs);
                ImGui::Checkbox("Use MS", &config.tSettings.useMs);
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Noise Settings"))
            {
                ImGui::InputInt("FFT Size", &config.nSettings.nfft);
                ImGui::InputInt("Hop Length", &config.nSettings.hopSize);
                ImGui::InputInt("Start sample", &config.nSettings.startSample);
                ImGui::Checkbox("Use previous lists", &config.nSettings.useList);
                ImGui::TreePop();
            }
        }

        // End the window
        ImGui::End();
    }
};