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
	NoiseTracker(const Sitrano::NoiseSettings& settings,
		const Sitrano::AnalysisUnit& ana, 
		const Sitrano::Results& r,
		float startFreq);
	~NoiseTracker();

	struct Band {
		double f_low;
		double f_high;
		double f_center;
	};
	void applyFrameTable(std::vector<int> table);
	std::vector<std::vector<float>> analyze();

private:

	const Sitrano::AnalysisUnit& unit;

	bool isInVector(float freq, const std::vector<Sitrano::Peak>& topFreqs);
	int findLargeBin(float freq, std::vector<Band> bands);
	float sr{ 96000.f };
	bool useList;
	float startFreq{ 20.f };
	const int N;
	const int num;
	const int hop;

	fftwf_plan plan;
	float* fft_in;
	fftwf_complex* fft_out;

	const std::vector<float>& sf;
	const std::vector<Sitrano::Peak>& topFreqs;
	std::vector<Band> bands;
	std::vector<int> frameTable;
	std::vector<std::vector<float>> results;
};