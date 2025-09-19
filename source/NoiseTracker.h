#pragma once

#include<iostream>
#include<vector>
#include<fftw3.h>
#include<cmath>
#include<sndfile.h>
#include<fstream>

class NoiseTracker {
public:
	NoiseTracker(float start, fftw_complex* fftOut) :
		out(fftOut),
		startFreq(start)
	{};
	~NoiseTracker() {};

	struct Band {
		double f_low;
		double f_high;
		double f_center;
	};



private:
	float startFreq{0.f};
	fftw_complex* out;

};