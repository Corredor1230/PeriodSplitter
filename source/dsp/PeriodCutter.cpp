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
    mStartSample(start),
    mSampleRate(unit.sampleRate)
{}

// This is the main public function that does all the work.
std::vector<uint32_t> PeriodCutter::findPeriodSamples()
{
    int startSample = findStartTransient();

    // dynamic values
    int expectedPeriodLength = static_cast<int>(mSampleRate / mPitch);
    int windowSize = (int)(((float)mSampleRate / mPitch) + 1);
    int expectedNumPeriods = static_cast<int>(mUnit.soundFile.size() / expectedPeriodLength);

    // Pre-find all zero crossings
    std::vector<int> zeroCrossings = findZeroCrossings(mUnit.soundFile, startSample);
    if (zeroCrossings.empty()) {
        std::cerr << "PeriodCutter: No zero crossings found. Aborting." << std::endl;
        return {};
    }

    // Find the initial period start
    int initialSample = startSample + static_cast<int>(mSampleRate * (mConfig.periodStartOffsetMs / 1000.0));
    int peakSample = findPeakSample(mUnit.soundFile, initialSample, initialSample + windowSize, true);
    int periodStart = findNearestZeroCached(zeroCrossings, peakSample);
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
                int nz = findNearestZeroCached(zeroCrossings, lastStart + expectedPeriodLength);
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
                int nz = findNearestZeroCached(zeroCrossings, filePos);
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
                int nz = findNearestZeroCached(zeroCrossings, filePos);
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


// --- Private Implementation Methods ---

int PeriodCutter::findStartTransient()
{
    std::vector<float> rmsList;
    bool transientFound = false;
    int tInitSample = 0;

    // Get parameters from config
    int rmsSize = (mConfig.transientRmsSizeMs / 1000.f) * mSampleRate;
    int rmsHopLength = std::max(1, (int)((float)rmsSize * mConfig.transientRmsHopRatio));
    float factor = mConfig.transientFactor;
    float threshold = mConfig.transientThreshold;

    std::vector<float> rmsWindow(rmsSize); // This will hold the transient

    for (int samp = mStartSample; samp < mUnit.soundFile.size(); samp += rmsHopLength)
    {
        float windowSum = 0.f;
        std::vector<float> tempRms(rmsSize);

        for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < mUnit.soundFile.size(); rmsSamp++)
        {
            float x = mUnit.soundFile[samp + rmsSamp];
            windowSum += (x * x);
            tempRms[rmsSamp] = x;
        }

        float rms = std::sqrt(windowSum / rmsSize);
        if (rmsList.size() < 2) {
            rmsList.push_back(rms);
            continue;
        }

        float rmsRatio = rmsList.back() == 0.f ? 1.f : rms / rmsList.back();
        rmsList.push_back(rms);

        if (rmsRatio > factor && rms > threshold) {
            transientFound = true;
            rmsWindow = tempRms; // Save the window that triggered
            tInitSample = samp;
            break;
        }
    }

    if (!transientFound) {
        // No transient, just find the first peak in a small window
        int peakSample = findPeakSample(mUnit.soundFile, 0, std::min((int)mUnit.soundFile.size(), 4096), true);
        return findPreviousZero(mUnit.soundFile, peakSample);
    }

    int rmsPeakSample = findPeak(rmsWindow);
    int peakSample = tInitSample + rmsPeakSample;
    return findPreviousZero(mUnit.soundFile, peakSample);
}


// --- Static Utility Functions (copied from your .cpp) ---
// (These are unchanged, just marked static)

int PeriodCutter::findPeak(const std::vector<float>& transientVector)
{
    if (transientVector.empty()) return 0;

    std::vector<float> absTransient(transientVector.size());
    for (int i = 0; i < absTransient.size(); i++)
    {
        absTransient[i] = std::abs(transientVector[i]);
    }

    int peakSample = std::distance(absTransient.begin(),
        std::max_element(absTransient.begin(), absTransient.end()));
    return peakSample;
}

int PeriodCutter::findPreviousZero(const std::vector<float>& signal, int startSample)
{
    int zeroSample = 0;
    for (int i = startSample; i > 2 && i < signal.size(); i--)
    {
        bool zeroFound = false;
        zeroFound = (signal[i] >= 0.0 && signal[i - 2] <= 0.0)
            || (signal[i] <= 0.0 && signal[i - 2] >= 0.0);

        if (zeroFound) {
            if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
                zeroSample = i - 2;
            else
                zeroSample = i - 1;
            break;
        }
    }
    return zeroSample;
}

int PeriodCutter::findNextZero(const std::vector<float>& signal, int startSample)
{
    int zeroSample = 0;
    for (int i = startSample; i > 1 && i < signal.size(); i++)
    {
        bool zeroFound = false;
        zeroFound = (signal[i] >= 0.0 && signal[i - 2] < 0.0)
            || (signal[i] <= 0.0 && signal[i - 2] > 0.0);

        if (zeroFound) {
            if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
                zeroSample = i - 2;
            else
                zeroSample = i - 1;
            break;
        }
    }
    return zeroSample;
}

int PeriodCutter::findNearestZero(const std::vector<float>& signal, int startSample)
{
    int prevZero = findPreviousZero(signal, startSample);
    int nextZero = findNextZero(signal, startSample);

    int distancePrev = std::abs(startSample - prevZero);
    int distanceNext = std::abs(startSample - nextZero);

    if (distancePrev < distanceNext)
        return prevZero;
    else
        return nextZero;
}

int PeriodCutter::findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute)
{
    int peakSample = startSample;
    float peakValue = -1.0e30f; // Use a very small number, not 0
    if (useAbsolute) peakValue = 0.f;

    int safeEndSample = std::min((int)signal.size(), endSample);
    int safeStartSample = std::max(0, startSample);

    for (int sample = safeStartSample; sample < safeEndSample; sample++)
    {
        if (useAbsolute) {
            if (std::abs(signal[sample]) > peakValue) {
                peakValue = std::abs(signal[sample]);
                peakSample = sample;
            }
        }
        else {
            if (signal[sample] > peakValue) {
                peakValue = signal[sample];
                peakSample = sample;
            }
        }
    }
    return peakSample;
}

std::vector<int> PeriodCutter::findZeroCrossings(const std::vector<float>& signal, int initSample)
{
    std::vector<int> zeroCrossings;
    for (int i = std::max(1, initSample); i < signal.size(); ++i)
    {
        if ((signal[i - 1] < 0 && signal[i] >= 0) ||
            (signal[i - 1] > 0 && signal[i] <= 0))
        {
            // Find the sample closer to zero
            zeroCrossings.push_back(std::abs(signal[i - 1]) < std::abs(signal[i]) ? (i - 1) : i);
        }
    }
    return zeroCrossings;
}

int PeriodCutter::findNearestZeroCached(const std::vector<int>& zeroCrossings, int sample)
{
    if (zeroCrossings.empty()) return sample; // Safety check

    auto it = std::lower_bound(zeroCrossings.begin(), zeroCrossings.end(), sample);
    if (it == zeroCrossings.end()) return zeroCrossings.back();
    if (it == zeroCrossings.begin()) return *it;

    int after = *it;
    int before = *(it - 1);
    return (std::abs(after - sample) < std::abs(before - sample)) ? after : before;
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
