#pragma once

#include"support/IWindowStrategy.h"

class AudioChunkStrategy : public IWindowStrategy {
public:


	LoopParameters getLoopParameters(const std::vector<uint32_t>& sampleList,
		const Sihat::AnalysisUnit& unit,
		const Sihat::AnalysisConfig& config, const int firstSample = 0) const override {
		int chunkSize = 0;
		if (firstSample > 0) chunkSize = unit.soundFile.size() - firstSample;
		else chunkSize = unit.soundFile.size() - sampleList[0];
		size_t numFrames = chunkSize / config.hopSize;
		return { numFrames, 1 };
	}

	bool processFrame(size_t i, 
		int frameStep, 
		const Sihat::AnalysisConfig& config,
		const Sihat::AnalysisUnit& unit,
		const std::vector<uint32_t>& sampleList,
		Sihat::HarmonicResults& results, 
		float* inputBuffer, 
		int& outPeriodLength,
		int firstSample) override 
	{
		int start = i * config.hopSize + firstSample;
		int end = start + config.nfft;
		if (end > unit.soundFile.size()) return false;

		results.finalSamples.push_back(start);
		for (int j = 0; j < config.nfft && start + j < unit.soundFile.size(); ++j) {
			inputBuffer[j] = unit.soundFile[start + j];
		}

		outPeriodLength = 0;
		return true;
	}
};