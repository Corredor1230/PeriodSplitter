#include "IWindowStrategy.h"

class PeriodLoopStrategy : public IWindowStrategy {
public:
    LoopParameters getLoopParameters(const std::vector<uint32_t>& sampleList,
        const Sitrano::AnalysisUnit& unit,
        const Sitrano::AnalysisConfig& config) const override {
        size_t numFrames = (sampleList.size() - 1) / 2;
        return { numFrames, 2 }; // {numFrames, frameStep}
    }

    bool processFrame(size_t i,
        int frameStep,
        const Sitrano::AnalysisConfig& config,
        const Sitrano::AnalysisUnit& unit,
        const std::vector<uint32_t>& sampleList,
        Sitrano::HarmonicResults& results,
        float* inputBuffer,
        int& outPeriodLength) override 

    {
        int start = sampleList[i];
        int end = sampleList[i + frameStep];
        outPeriodLength = end - start;

        if (outPeriodLength <= 0 || end > unit.soundFile.size()) return false;

        for (int j = 0; j < config.nfft; ++j) {
            int srcIdx = start + (j % outPeriodLength);
            inputBuffer[j] = (srcIdx < unit.soundFile.size()) ? unit.soundFile[srcIdx] : 0.f;
        }

        // This logic is also encapsulated here
        for (int b = 0; b < frameStep && i + b < sampleList.size(); b++) {
            results.finalSamples.push_back(sampleList[i + b]);
        }
        return true;
    }
};