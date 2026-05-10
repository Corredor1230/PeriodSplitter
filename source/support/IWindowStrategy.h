#pragma once

#include"include/SitranoHeader.h"

struct LoopParameters {
	size_t numFrames;
	int frameStep;
};

class IWindowStrategy {
public:
	virtual ~IWindowStrategy() = default;
	virtual LoopParameters getLoopParameters(const std::vector<uint32_t>& sampleList,
		const Sihat::AnalysisUnit& unit, 
		const Sihat::AnalysisConfig& config, const int firstSample) const = 0;
	virtual bool processFrame(size_t i,
		int frameStep,
		const Sihat::AnalysisConfig& config,
		const Sihat::AnalysisUnit& unit,
		const std::vector<uint32_t>& sampleList, 
		Sihat::HarmonicResults& results,
		float* inputBuffer,
		int& outPeriodLength,
		int firstSample) = 0;
};