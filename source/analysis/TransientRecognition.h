#pragma once

#include<vector>
#include"include/SitranoHeader.h"
#include<fftw3.h>

class Transient
{
public:


	Transient(
		const Sihat::AnalysisUnit& unit,
		const Sihat::TransientSettings& conf,
		const Sihat::TransientFFTSettings& ffts,
		const float pitch
	);

	~Transient();

	/**
	 * @brief Finds the first major transient in the signal based on RMS.
	 * @return The sample index of the nearest zero-crossing *before* the transient.
	 */
	Sihat::TransientResults findStartTransient();

private:
	/*
	* These are used to find the RMS transient
	*/
	const float pitch;
	const int correlationOffset;
	const float corrThreshold;

	int rmsSize;
	int rmsHopLength;
	int preAttack;
	float factor;
	float threshold;
	float sampleRate;
	const Sihat::TransientSettings& tSettings;
	const std::vector<float>& aud;
	Sihat::SampleRange findFromRMS();
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
	const int waveletNumber = 64;
	const float flatnessThresh;
	const bool useFFT;
	Sihat::SampleRange findFromFFT();
	Sihat::SampleRange findWithCrossCorrelation(int offset, int firstSample);
};