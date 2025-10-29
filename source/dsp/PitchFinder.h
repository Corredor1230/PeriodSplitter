#pragma once

#include <memory> // For std::unique_ptr
#include "include/SitranoHeader.h" // For Sitrano types

class PyinCpp;

class PitchFinder {
public:
    PitchFinder(const Sitrano::AnalysisUnit& u,
        const Sitrano::AnalysisConfig& config);

    ~PitchFinder();

    float findPitch();

private:
    /**
     * @brief Finds the densest cluster in a 1D dataset and returns its average.
     * @param data The vector of pitches. Will be sorted in-place (moved).
     * @param threshold The max distance between values in a cluster.
     * @param minFreq Values below this are ignored.
     * @param maxFreq Values above this are ignored.
     * @return The average of the densest cluster, or NAN if no valid cluster found.
     */
    static float findMode(std::vector<float> data, float threshold,
        float minFreq, float maxFreq);

    // --- Member Variables ---

    const Sitrano::AnalysisUnit& unit;
    const Sitrano::AnalysisConfig& mConfig; // Store the config
    float tolerance = 10.f;

    std::unique_ptr<PyinCpp> mPitchDetector;
};
