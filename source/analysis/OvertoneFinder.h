#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include"include/SitranoHeader.h"

class OvertoneFinder {
public:
	OvertoneFinder(const Sitrano::AnalysisUnit& unit,
		const Sitrano::AnalysisConfig& config);
	~OvertoneFinder() {
		fftwf_destroy_plan(plan);
	}

	std::vector<Sitrano::Peak> getRelevantOvertones(const std::vector<double>& checkSignal,
		float pitch);

private:
	const Sitrano::AnalysisUnit& unit;
	const Sitrano::OvertoneSettings& settings;
	const Sitrano::AnalysisConfig& config;
	void initFFTW();
	float* input;
	int N;
	int Nout;
	bool absThreshold;
	float threshold;
	fftwf_plan plan;
	fftwf_complex* output;
};