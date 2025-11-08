#include "PitchFinder.h"

// --- Implementation-Only Includes ---
// 1. All "heavy" or "detail" headers are now hidden in the .cpp.
// Anyone including PitchFinder.h doesn't pay this compile cost.
#include <libpyincpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath> // For NAN

// 2. The constructor initializes the Pimpl and stores the config
PitchFinder::PitchFinder(const Sitrano::AnalysisUnit& u,
    const Sitrano::AnalysisConfig& config) :
    unit(u),
    mConfig(config),
    mPitchDetector(std::make_unique<PyinCpp>(u.sampleRate, mConfig.nfft, mConfig.nfft / 2))
{
    // The mPitchDetector is now created *once* and stored.
    tolerance = mConfig.tolerance;
}

// 3. The destructor is defined here, where PyinCpp is a complete type.
// This is required for std::unique_ptr to work with the forward declaration.
PitchFinder::~PitchFinder() = default;


float PitchFinder::findPitch()
{
    // 4. Use the *member* pitch detector. No more re-creating it.
    std::vector<float> pitches = mPitchDetector->feed(unit.soundFile);

    // 5. Get parameters from the config, not hard-coded magic numbers.
    float foundPitch = findMode(std::move(pitches), // 6. Move pitches (avoids copy)
        mConfig.pSettings.modeThreshold,
        mConfig.pSettings.minFreq,
        mConfig.pSettings.maxFreq);

    if (std::isnan(foundPitch)) {
        std::cerr << "PitchFinder: Could not find a valid pitch mode." << std::endl;
        // Fallback to metadata pitch if analysis fails
        return Sitrano::getPitchFromFilename(unit.filename);
    }

    float metaPitch = Sitrano::getPitchFromFilename(unit.filename);
    float pitch{ 0.0 };

    if (metaPitch <= 0.0) {
        pitch = foundPitch; // No metadata, trust the analysis
    }
    else {
        // Use the config for tolerance
        float tolerance = Sitrano::cents_to_hz(metaPitch, mConfig.pSettings.toleranceInCents);

        // Adjudicate between analyzed pitch and metadata pitch
        if (foundPitch > metaPitch - tolerance && foundPitch < metaPitch + tolerance)
            pitch = foundPitch; // Analysis matches metadata, great!
        else
            pitch = metaPitch; // Mismatch, trust the metadata
    }

    // You can control this with a verbose flag in mConfig
    if (mConfig.verbose) {
        std::cout << "PitchFinder: Analyzed pitch = " << foundPitch << " Hz" << std::endl;
        std::cout << "PitchFinder: Metadata pitch = " << metaPitch << " Hz" << std::endl;
        std::cout << "PitchFinder: Final pitch = " << pitch << " Hz\n";
    }

    return pitch;
}

// 7. 'findMode' is now static as it doesn't use any member data.
// 8. It takes 'data' by VALUE, allowing us to sort it in-place.
float PitchFinder::findMode(std::vector<float> data, float threshold,
    float minFreq, float maxFreq)
{
    if (data.empty()) return NAN;

    // Step 1: Sort the data (this is our *own copy* thanks to pass-by-value)
    std::sort(data.begin(), data.end());

    // Step 2: Slide a window
    int maxCount = 0;
    float modeSum = 0.0f;
    int modeCount = 0;

    for (size_t i = 0; i < data.size(); ++i) {
        float ref = data[i];

        if (ref < minFreq || ref > maxFreq)
            continue;

        float sum = 0.0f;
        int count = 0;

        // Include all values within [ref, ref + threshold)
        for (size_t j = i; j < data.size() && data[j] - ref < threshold; ++j) {
            sum += data[j];
            ++count;
        }

        if (count > maxCount) {
            maxCount = count;
            modeSum = sum;
            modeCount = count;
        }

        // Small optimization: if a cluster has maxCount,
        // no cluster starting *inside* it can be larger.
        // We can skip ahead to the end of this cluster.
        if (count > 0) {
            i += count - 1;
        }
    }

    return modeCount > 0 ? (modeSum / modeCount) : NAN;
}