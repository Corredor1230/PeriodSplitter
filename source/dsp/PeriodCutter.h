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
        float pitch, int start);

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
    const int mStartSample;
    const float mPitch;
    const float mSampleRate;

    // --- Private Implementation Methods ---

    /**
     * @brief Finds the first major transient in the signal based on RMS.
     * @return The sample index of the nearest zero-crossing *before* the transient.
     */
    int findStartTransient();

    // --- Static Utility Functions (Helpers) ---
    // These are static because they don't depend on the object's state,
    // only on their arguments.

    static int findPeak(const std::vector<float>& transientVector);

    static int findPreviousZero(const std::vector<float>& signal, int startSample);

    static int findNextZero(const std::vector<float>& signal, int startSample);

    static int findNearestZero(const std::vector<float>& signal, int startSample);

    static int findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute = true);

    static std::vector<int> findZeroCrossings(const std::vector<float>& signal, int initSample);

    static int findNearestZeroCached(const std::vector<int>& zeroCrossings, int sample);

    static float signalCorrelation(const std::vector<float>& window, const std::vector<float>& signal, int startSample);

    struct CorrelationState {
        float numerator = 0.f;
        float squareB = 0.f;
    };

    static float signalCorrelationRolling(const std::vector<float>& window, float squareA,
        const std::vector<float>& signal, int startSample,
        CorrelationState& state);
};
