#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include"include/SitranoHeader.h"

class OvertoneFinder {
public:
	OvertoneFinder(const Sihat::AnalysisUnit& unit,
		const Sihat::AnalysisConfig& config);
	~OvertoneFinder() {
		fftwf_destroy_plan(plan);
	}

	std::vector<Sihat::Peak> getRelevantOvertones(const std::vector<double>& checkSignal,
		float pitch);

private:
	const Sihat::AnalysisUnit& unit;
	const Sihat::OvertoneSettings& settings;
	const Sihat::AnalysisConfig& config;
	void initFFTW();
	float* input;
	int N;
	int Nout;
	bool absThreshold;
	float threshold;
	fftwf_plan plan;
	fftwf_complex* output;
};