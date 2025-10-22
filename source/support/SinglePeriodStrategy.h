#include "IWindowStrategy.h"

class SinglePeriodStrategy : public IWindowStrategy {
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
            int srcIdx = start + j;

        return true;
    }
};