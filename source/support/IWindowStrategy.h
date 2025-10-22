#pragma once

#include"include/SitranoHeader.h"

struct LoopParameters {
	size_t numFrames;
	int frameStep;
};

class IWindowStrategy {
public:
	virtual ~IWindowStrategy() = default;
	virtual LoopParameters getLoopParameters(const Sitrano::Results& results,
		const Sitrano::AnalysisUnit& unit) const = 0;
	virtual bool processFrame(size_t i,
		int frameStep,
		const Sitrano::AnalysisUnit& unit,
		Sitrano::Results& results,
		float* inputBuffer,
		int& outPeriodLength) = 0;
};