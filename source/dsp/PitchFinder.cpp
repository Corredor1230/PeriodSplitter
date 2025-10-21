#include"PitchFinder.h"

PitchFinder::PitchFinder(const Sitrano::AnalysisUnit& u) :
	unit(u)
{

}

float PitchFinder::getPitchFromFilename(const std::string& filename)
{
    std::string name = filename;
    std::vector<std::string> parts;
    size_t start = 0, pos;

    while ((pos = name.find('_', start)) != std::string::npos) {
        parts.push_back(name.substr(start, pos - start));
        start = pos + 1;
    }
    parts.push_back(name.substr(start));

    float outPitch = std::stof(parts[2]);

    return outPitch;
}

float PitchFinder::findPitch()
{
    PyinCpp pitchDetector(unit.sampleRate);

    std::vector<float> pitches = pitchDetector.feed(unit.soundFile);

    float foundPitch = findMode(pitches);

    //std::cout << "Current pitch: " << pitch << "Hz\n";

    float metaPitch = getPitchFromFilename(unit.filename);
    float pitch{ 0.0 };
    float tolerance = 3.0;

    if (foundPitch > metaPitch - tolerance && foundPitch < metaPitch + tolerance)
        pitch = foundPitch;
    else
        pitch = metaPitch;

    std::cout << "Current pitch: " << pitch << 'Hz\n';

    return pitch;
}

float PitchFinder::findMode(const std::vector<float>& data, float threshold,
	float minFreq, float maxFreq)
{
    if (data.empty()) return NAN;

    // Step 1: Sort the data
    std::vector<float> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Step 2: Slide a window and track the group with the highest count
    int maxCount = 0;
    float modeSum = 0.0f;
    int modeCount = 0;

    for (size_t i = 0; i < sortedData.size(); ++i) {
        float ref = sortedData[i];
        float sum = 0.0f;
        int count = 0;

        if (ref < minFreq || ref > maxFreq)
            continue;

        // Include all values within [ref, ref + threshold)
        for (size_t j = i; j < sortedData.size() && sortedData[j] - ref < threshold; ++j) {
            sum += sortedData[j];
            ++count;
        }

        if (count > maxCount) {
            maxCount = count;
            modeSum = sum;
            modeCount = count;
        }
    }

    float mode = modeSum / modeCount;

    return modeCount > 0 ? mode : NAN;
}