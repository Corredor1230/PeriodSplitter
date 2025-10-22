#pragma once

#include"support/IWindowStrategy.h"

class AudioChunkStrategy : public IWindowStrategy {
public:
	LoopParameters getLoopParameters(const Sitrano::Results& results,
		const Sitrano::AnalysisUnit& unit) const override {
		int chunkSize = unit.soundFile.size() - results.sampleList[0];
		size_t numFrames = chunkSize / unit.hopSize;
		return { numFrames, 1 };
	}

	bool processFrame(size_t i, int frameStep, const Sitrano::AnalysisUnit& unit,
		Sitrano::Results& results, float* inputBuffer, int& outPeriodLength) override {
		int start = i * unit.hopSize;
		int end = start + unit.nfft;
		if (end > unit.soundFile.size()) return false;

		results.finalSamples.push_back(start);
		for (int j = 0; j < unit.nfft; ++j) {
			inputBuffer[j] = unit.soundFile[start + j];
		}

		outPeriodLength = 0;
		return true;
	}
};