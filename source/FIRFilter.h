#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

constexpr double PI = 3.14159265358979323846;

class FIRHighPass {
public:
    FIRHighPass(float sampleRate, float cutoffFreq, float transitionWidth);

    // Process a single sample
    double processSample(float input);

    // Process a block of samples
    std::vector<float> processBuffer(const std::vector<float>& input);
    void processAudioFile(std::vector<float>& input, const int bufferSize);

    // Access the filter kernel (for visualization or testing)
    const std::vector<float>& getTaps() const;

private:
    std::vector<float> taps;
    std::vector<float> history;

    void designFilter(float sampleRate, float cutoffFreq, float transitionWidth);
};
