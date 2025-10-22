#include "IWindowStrategy.h"

class PeriodLoopStrategy : public IWindowStrategy {
public:
    LoopParameters getLoopParameters(const Sitrano::Results& results, const Sitrano::AnalysisUnit& unit) const override {
        size_t numFrames = (results.sampleList.size() - 1) / 2;
        return { numFrames, 2 }; // {numFrames, frameStep}
    }

    bool processFrame(size_t i, int frameStep, const Sitrano::AnalysisUnit& unit, Sitrano::Results& results,
        float* inputBuffer, int& outPeriodLength) override {

        int start = results.sampleList[i];
        int end = results.sampleList[i + frameStep];
        outPeriodLength = end - start;

        if (outPeriodLength <= 0 || end > unit.soundFile.size()) return false;

        for (int j = 0; j < unit.nfft; ++j) {
            int srcIdx = start + (j % outPeriodLength);
            inputBuffer[j] = (srcIdx < unit.soundFile.size()) ? unit.soundFile[srcIdx] : 0.f;
        }

        // This logic is also encapsulated here
        for (int b = 0; b < frameStep && i + b < results.sampleList.size(); b++) {
            results.finalSamples.push_back(results.sampleList[i + b]);
        }
        return true;
    }
};