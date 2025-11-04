#pragma once

#include <vector>
#include <cstdint> // For uint32_t
#include "include/SitranoHeader.h"

class PeriodCutter
{
public:
    /**
     * @brief Constructs a stateless PeriodCutter.
     * @param unit The input audio data.
     * @param config The configuration for transient and correlation.
     * @param pitch The fundamental pitch (from PitchFinder) to search for.
     */
    PeriodCutter(const Sitrano::AnalysisUnit& unit,
        const Sitrano::CorrelationSettings& config,
        float pitch, const Sitrano::SampleRange& transient);

    ~PeriodCutter() = default;

    /**
     * @brief Performs the full analysis to find period zero-crossings.
     * This is the main entry point for the class.
     * @return A vector of sample indices (uint32_t) for each period start.
     */
    std::vector<uint32_t> findPeriodSamples();

private:
    // --- Private Member Data ---
    const Sitrano::AnalysisUnit& mUnit;
    const Sitrano::CorrelationSettings& mConfig;
    const int startSample;
    const float mPitch;
    const float mSampleRate;
    const int tStart;

    // --- Private Implementation Methods ---

    static float signalCorrelation(const std::vector<float>& window, const std::vector<float>& signal, int startSample);

    struct CorrelationState {
        float numerator = 0.f;
        float squareB = 0.f;
    };

    static float signalCorrelationRolling(const std::vector<float>& window, float squareA,
        const std::vector<float>& signal, int startSample,
        CorrelationState& state);
};
