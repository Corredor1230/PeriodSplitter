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
	NoiseTracker(const Sihat::NoiseSettings& settings,
		const Sihat::AnalysisUnit& ana, 
		const Sihat::Results& r,
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

	const Sihat::AnalysisUnit& unit;

	bool isInTop(float freq, const std::vector<Sihat::Peak>& topFreqs,
		const std::vector<Band>& bands);
	int findLargeBin(float freq, std::vector<Band> bands);
	float sr{ 96000.f };
	bool useList;
	float startFreq{ 20.f };
	const int N;
	const int num;
	const float minFreq;
	const float maxFreq;
	const int hop;

	fftwf_plan plan;
	float* fft_in;
	fftwf_complex* fft_out;

	const std::vector<float>& sf;
	const std::vector<Sihat::Peak>& topFreqs;
	std::vector<Band> bands;
	std::vector<int> frameTable;
	std::vector<std::vector<float>> results;
};