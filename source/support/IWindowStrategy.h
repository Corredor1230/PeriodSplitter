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
		const Sitrano::AnalysisUnit& unit, 
		const Sitrano::AnalysisConfig& config) const = 0;
	virtual bool processFrame(size_t i,
		int frameStep,
		const Sitrano::AnalysisConfig& config,
		const Sitrano::AnalysisUnit& unit,
		const std::vector<uint32_t>& sampleList, 
		Sitrano::HarmonicResults& results,
		float* inputBuffer,
		int& outPeriodLength) = 0;
};