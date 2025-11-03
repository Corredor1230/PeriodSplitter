#include "PeriodCutter.h"

#include <iostream>
#include <deque>
#include <array>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <future>   // For std::async
#include <chrono>   // For std::chrono
#include <numeric>  // For std::distance

PeriodCutter::PeriodCutter(const Sitrano::AnalysisUnit& unit,
    const Sitrano::CorrelationSettings& config,
    float pitch, int start) :
    mUnit(unit),
    mConfig(config),
    mPitch(pitch),
    startSample(start),
    mSampleRate(unit.sampleRate)
{}

// This is the main public function that does all the work.
std::vector<uint32_t> PeriodCutter::findPeriodSamples()
{

    // dynamic values
    int expectedPeriodLength = static_cast<int>(mSampleRate / mPitch);
    int windowSize = (int)(((float)mSampleRate / mPitch) + 1);
    int expectedNumPeriods = static_cast<int>(mUnit.soundFile.size() / expectedPeriodLength);

    // Pre-find all zero crossings
    std::vector<int> zeroCrossings = Sitrano::findZeroCrossings(mUnit.soundFile, startSample);
    if (zeroCrossings.empty()) {
        std::cerr << "PeriodCutter: No zero crossings found. Aborting." << std::endl;
        return {};
    }

    // Find the initial period start
    int initialSample = startSample + static_cast<int>(mSampleRate * (mConfig.periodStartOffsetMs / 1000.0));
    int peakSample = Sitrano::findPeakSample(mUnit.soundFile, initialSample, initialSample + windowSize, true);
    int periodStart = Sitrano::findNearestCachedZero(zeroCrossings, peakSample);
    int margin = std::max<int>(4, std::ceil<int>((float)expectedPeriodLength * 0.02));

    // create the window template
    std::vector<float> window(windowSize);
    for (int i = 0; i < windowSize; i++)
        window[i] = (periodStart + i < mUnit.soundFile.size()) ? mUnit.soundFile[periodStart + i] : 0.f;
    std::vector<float> backWindow = window;

    // Precompute window norm
    float squareA = 0.f;
    for (float w : window) squareA += w * w;
    if (squareA == 0.f) {
        std::cerr << "PeriodCutter: Initial window is silent. Aborting." << std::endl;
        return {}; // Avoid divide by zero
    }

    std::deque<int> forwardList;
    std::deque<int> backwardList;
    forwardList.push_back(periodStart);
    backwardList.push_back(periodStart);

    // --- Define the parallel search lambdas ---
    auto forwardSearch = [&]() {
        std::vector<float> corrVals;
        CorrelationState state;
        int filePos = periodStart;
        bool first = true;
        int lastStart = filePos;

        while (filePos < mUnit.soundFile.size()) {
            if (first) {
                lastStart = filePos;
                // Re-calculate window and normA for the new starting point
                squareA = 0.f;
                for (int i = 0; i < windowSize; i++) {
                    window[i] = (lastStart + i < mUnit.soundFile.size()) ? mUnit.soundFile[lastStart + i] : 0.f;
                    squareA += window[i] * window[i];
                }
                if (squareA == 0.f) squareA = 0.0001f; // Avoid div by zero

                first = false;
                filePos += (expectedPeriodLength - margin - 2);
            }

            float corr = signalCorrelationRolling(window, squareA, mUnit.soundFile, filePos, state);
            if (corr > mConfig.correlationThreshold) {
                corrVals.push_back(corr);
            }
            else if (filePos > lastStart + expectedPeriodLength + 3) {
                int nz = Sitrano::findNearestCachedZero(zeroCrossings, lastStart + expectedPeriodLength);
                if (nz == forwardList.back()) { filePos++; continue; }
                if (nz - forwardList.back() >= expectedPeriodLength - margin &&
                    nz - forwardList.back() <= expectedPeriodLength + margin)
                {
                    forwardList.push_back(nz);
                    filePos = nz;
                }
                else
                {
                    forwardList.push_back(lastStart + expectedPeriodLength);
                    filePos = lastStart + expectedPeriodLength;
                }
                first = true;
            }
            else {
                filePos++;
                continue;
            }

            if (corrVals.size() < 3) continue;
            int lc = corrVals.size() - 1;
            bool isPeak = (corrVals[lc - 1] > corrVals[lc]) && (corrVals[lc - 1] > corrVals[lc - 2]);

            if (isPeak) {
                int nz = Sitrano::findNearestCachedZero(zeroCrossings, filePos);
                if (nz == forwardList.back()) { filePos++; continue; }
                int currentWindowSize = nz - forwardList.back();
                if (currentWindowSize >= (expectedPeriodLength - margin))
                {
                    forwardList.push_back(nz);
                    filePos = nz;
                }
                else
                {
                    int temporaryNewStart = forwardList.back() + expectedPeriodLength;
                    forwardList.push_back(temporaryNewStart);
                    filePos = temporaryNewStart;
                }
                first = true;
            }
            if (forwardList.size() > expectedNumPeriods * 1.1) break;
            filePos++;
        }
        };

    auto backwardSearch = [&]() {
        std::vector<float> corrVals;
        int filePos = periodStart;
        bool first = true;
        int lastStart = filePos;
        float backSquareA = squareA; // Use the initial squareA

        while (filePos > (startSample - windowSize / 2)) {
            if (first) {
                first = false;
                lastStart = filePos;
                // Recalculate backWindow and its norm
                backSquareA = 0.f;
                for (int i = 0; i < windowSize; i++) {
                    backWindow[i] = (lastStart + i) > 1 ? mUnit.soundFile[lastStart + i] : 0.f;
                    backSquareA += backWindow[i] * backWindow[i];
                }
                if (backSquareA == 0.f) backSquareA = 0.0001f;

                filePos -= (expectedPeriodLength - margin);
                corrVals.push_back(0.0);
            }

            float corr = signalCorrelation(backWindow, mUnit.soundFile, filePos);
            corrVals.push_back(corr);

            if (corrVals.size() < 3) {
                filePos--;
                continue;
            }
            int lc = corrVals.size() - 1;

            bool aboveThreshold = (corrVals[lc] > mConfig.correlationThreshold);
            bool isPeak = (corrVals[lc - 1] > corrVals[lc]) && (corrVals[lc - 1] > corrVals[lc - 2]);

            if (isPeak && aboveThreshold) {
                int nz = Sitrano::findNearestCachedZero(zeroCrossings, filePos);
                int currentWindowSize = backwardList.front() - nz;
                if (currentWindowSize <= expectedPeriodLength + margin &&
                    currentWindowSize >= expectedPeriodLength - margin)
                {
                    backwardList.push_front(nz);
                    filePos = nz + 1;
                }
                else
                {
                    backwardList.push_front(lastStart - expectedPeriodLength);
                    filePos = lastStart - expectedPeriodLength;
                }
                first = true;
            }
            filePos--;
        }
        };

    // Run searches in parallel
    auto forwardFuture = std::async(std::launch::async, forwardSearch);
    auto backwardFuture = std::async(std::launch::async, backwardSearch);
    forwardFuture.get();
    backwardFuture.get();

    // Merge and return results
    std::vector<uint32_t> periodZeroes;
    periodZeroes.insert(periodZeroes.end(), backwardList.begin(), backwardList.end());
    periodZeroes.insert(periodZeroes.end(), ++forwardList.begin(), forwardList.end());

    return periodZeroes;
}

float PeriodCutter::signalCorrelation(const std::vector<float>& window, const std::vector<float>& signal, int startSample)
{
    float squareA = 0.f;
    float squareB = 0.f;
    float numerator = 0.f;

    for (int i = 0; i < window.size(); i++)
    {
        float a = window[i];
        float b = (i + startSample < signal.size() && i + startSample >= 0) ? signal[i + startSample] : 0.f;

        numerator += a * b;
        squareA += a * a;
        squareB += b * b;
    }
    float denominator = std::sqrt(squareA * squareB);
    return (denominator != 0.f) ? numerator / denominator : 0.f;
}

float PeriodCutter::signalCorrelationRolling(const std::vector<float>& window, float squareA,
    const std::vector<float>& signal, int startSample, CorrelationState& state)
{
    // Note: The 'first' flag was removed, this function is assumed to be
    // called in a rolling fashion, but the old one re-calculated everything anyway.
    // This is the same as the non-rolling version.
    state.numerator = 0.f;
    state.squareB = 0.f;

    for (int i = 0; i < window.size(); ++i)
    {
        float a = window[i];
        float b = (i + startSample < signal.size() && i + startSample >= 0) ? signal[i + startSample] : 0.f;
        state.numerator += a * b;
        state.squareB += b * b;
    }

    float denominator = std::sqrt(squareA * state.squareB);
    return (denominator != 0.f) ? state.numerator / denominator : 0.f;
}
