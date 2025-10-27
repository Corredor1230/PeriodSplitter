#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include<cmath>
#include<sndfile.h>
#include<fstream>
#include"include/SitranoHeader.h"

class NoiseTracker {
public:
	NoiseTracker(const Sitrano::AnalysisConfig& config,
		const Sitrano::AnalysisUnit& ana, 
		Sitrano::Results& r,
		float startFreq);
	~NoiseTracker();

	struct Band {
		double f_low;
		double f_high;
		double f_center;
	};
	void applyFrameTable(std::vector<int> table);
	void analyze();

private:

	const Sitrano::AnalysisUnit& unit;
	const Sitrano::AnalysisConfig& config;
	Sitrano::Results& r;

	bool useFrameTable = false;

	float startFreq{0.f};
	float sr{ 96000.f };
	int N;

	fftwf_plan plan;
	float* fft_in;
	fftwf_complex* fft_out;

	std::vector<float> sf;
	std::vector<Band> bands;
	std::vector<int> frameTable;
	std::vector<std::vector<float>> results;
};