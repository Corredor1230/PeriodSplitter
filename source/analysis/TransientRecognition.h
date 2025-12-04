#pragma once

#include<vector>
#include"include/SitranoHeader.h"
#include<fftw3.h>

class Transient
{
public:


	Transient(
		const Sitrano::AnalysisUnit& unit,
		const Sitrano::TransientSettings& conf,
		const Sitrano::TransientFFTSettings& ffts,
		const float pitch
	);

	/**
	 * @brief Finds the first major transient in the signal based on RMS.
	 * @return The sample index of the nearest zero-crossing *before* the transient.
	 */
	Sitrano::SampleRange findStartTransient();

private:
	/*
	* These are used to find the RMS transient
	*/
	const float pitch;
	const int correlationOffset;
	int rmsSize;
	int rmsHopLength;
	int preAttack;
	float factor;
	float threshold;
	float sampleRate;
	const Sitrano::TransientSettings& tSettings;
	const std::vector<float>& aud;
	Sitrano::SampleRange findFromRMS();
	int findFirstAboveThreshold(int startSample, float thresh);


	/*
	* These are used to find the FFT transient
	*/
	fftwf_plan plan;
	float* input;
	fftwf_complex* output;
	void initFFTW();
	const int nfft;
	const float sr;
	const int hop;
	const float flatnessThresh;
	const bool useFFT;
	Sitrano::SampleRange findFromFFT();
	Sitrano::SampleRange findWithCrossCorrelation(int offset, int firstSample);
};