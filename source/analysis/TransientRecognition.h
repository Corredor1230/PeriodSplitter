#pragma once

#include<vector>
#include"include/SitranoHeader.h"

class Transient
{
public:
	Transient(
		const Sitrano::AnalysisUnit& unit,
		const Sitrano::TransientSettings& conf
	);

	/**
	 * @brief Finds the first major transient in the signal based on RMS.
	 * @return The sample index of the nearest zero-crossing *before* the transient.
	 */
	Sitrano::SampleRange findStartTransient();

private:
	int rmsSize;
	int rmsHopLength;
	int backwards;
	float factor;
	float threshold;
	float sampleRate;
	const Sitrano::TransientSettings& tSettings;
	const std::vector<float>& aud;
};