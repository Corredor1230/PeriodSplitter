#pragma once

#include"support/IWindowStrategy.h"

class AudioChunkStrategy : public IWindowStrategy {
public:
	LoopParameters getLoopParameters(const std::vector<uint32_t>& sampleList,
		const Sitrano::AnalysisUnit& unit,
		const Sitrano::AnalysisConfig& config) const override {
		int chunkSize = unit.soundFile.size() - sampleList[0];
		size_t numFrames = chunkSize / config.hopSize;
		return { numFrames, 1 };
	}

	bool processFrame(size_t i, 
		int frameStep, 
		const Sitrano::AnalysisConfig& config,
		const Sitrano::AnalysisUnit& unit,
		const std::vector<uint32_t>& sampleList,
		Sitrano::HarmonicResults& results, 
		float* inputBuffer, 
		int& outPeriodLength) override {
		int start = i * config.hopSize;
		int end = start + config.nfft;
		if (end > unit.soundFile.size()) return false;

		results.finalSamples.push_back(start);
		for (int j = 0; j < config.nfft; ++j) {
			inputBuffer[j] = unit.soundFile[start + j];
		}

		outPeriodLength = 0;
		return true;
	}
};