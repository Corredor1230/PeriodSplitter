#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include<cmath>
#include<sndfile.h>
#include<fstream>

class NoiseTracker {
public:
	NoiseTracker(float start, float sr);
	~NoiseTracker() {};

	struct Band {
		double f_low;
		double f_high;
		double f_center;
	};
	void applyFrameTable(std::vector<int> table);
	void analyze(fftw_complex* fftOut) {};

private:
	bool useFrameTable = false;

	float startFreq{0.f};
	float sr{ 96000.f };

	std::vector<Band> bands;
	std::vector<int> frameTable;
	std::vector<std::vector<float>> results;
};