#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include<cmath>
#include<sndfile.h>
#include<fstream>

class NoiseTracker {
public:
	NoiseTracker(std::vector<float> audio, float start, float sr, int fftSize, int hopSize);
	~NoiseTracker() {};

	struct Band {
		double f_low;
		double f_high;
		double f_center;
	};
	void applyFrameTable(std::vector<int> table);
	void analyze();

private:
	int N{ 0 };
	int hop{ 0 };

	bool useFrameTable = false;

	float startFreq{0.f};
	float sr{ 96000.f };

	std::vector<float> sf;
	std::vector<Band> bands;
	std::vector<int> frameTable;
	std::vector<std::vector<float>> results;
};