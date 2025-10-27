#pragma once

#include<libpyincpp.h>
#include<iostream>
#include<string.h>
#include<vector>
#include<algorithm>
#include"include/SitranoHeader.h"

class PitchFinder {
public:

	PitchFinder(const Sitrano::AnalysisUnit& u);
	~PitchFinder() {};

	float findPitch();

private:

	float findMode(const std::vector<float>& data, float threshold = 1.0f,
		float minFreq = 60.0, float maxFreq = 1200.0);


	const Sitrano::AnalysisUnit& unit;
	std::string divChar;

};