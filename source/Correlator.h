#pragma once

#include<fftw3.h>
#include<libpyincpp.h>
#include<sndfile.h>
#include<iostream>
#include<vector>
#include<deque>
#include<array>
#include<string>
#include<fstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<future>
#include<chrono>

class Correlator
{
public:
	Correlator(std::vector<float>& file, SF_INFO& info, float sizeInMs = 50.f, float hopInMs = 0.25);
	~Correlator() {};

	void initialize(float pitch, float threshold = 0.1);
	void setWindowSize(int sizeInSamples = 512);
	void setWindowSize(float sizeInMs);
	void setThreshold(float thresh);
	void setPitch(float pitch);
	void enableHopping(bool hopping);
	void setHopLength(float hopInMs);
	void setHopLength(int hopInSamples);
	void printZeroList(std::string& filename);
	float findPitch(const std::vector<float>& signal);
	std::vector<int> getCorrelationPeaks();
	std::vector<int> getCorrelationZeroes();

	std::vector<std::vector<int>>	chopPositions;
	std::vector<std::vector<float>>	multiChannel;

private:

	struct CorrelationState
	{
		float numerator = 0.f;
		float squareB = 0.f;
	};

	std::vector<float>&	audioFile;
	std::vector<int>	peakList;
	std::vector<int>	zeroList;
	std::vector<int>	allZeroes;
	SF_INFO&			sfInfo;
	PyinCpp				pitchDetector;

	int		findStartTransient(int startSample, std::vector<float>& vector, 
			int rmsSize, int rmsHopLength, float transFactor = 2.0, float transThresh = 0.2);
	int		findStartTransient(int startSample, std::vector<float>& vector, 
			float rmsSizeInMs, float rmsHopLengthInRatio, float transFactor = 2.0, float transThresh = 0.2);

	int		findPeak(int startSample, std::vector<float> vector);

	int		findPreviousZero(std::vector<float>& signal, int startSample);
	int		findNextZero(std::vector<float>& signal, int startSample);
	int		findNearestZero(std::vector<float>& signal, int startSample);
	float	findPeakValue(std::vector<float>& signal, bool useAbsolute = true);
	int		findPeakSample(std::vector<float>& signal, int startSample, int endSample, bool useAbsolute = true);
	int		findPeakSample(std::vector<float>& signal, bool useAbsolute = true);
	std::vector<int> findPeriodSamples(std::vector<float>& signal, int startSample, 
				float msOffset, float inPitch, float correlationThreshold = 0.4);
	std::vector<int> findZeroCrossings(const std::vector<float>& signal, int initSample);
	int findNearestZeroCached(const std::vector<int>& zeroCrossings, int sample);
	float	findMode(const std::vector<float>& vctr, float threshold, float minFreq, float maxFreq);
	float	signalCorrelation(std::vector<float>& window, std::vector<float>& signal, int startSample);
	float	signalCorrelationRolling(const std::vector<float>& window, float squareA,
		const std::vector<float>& signal, int startSample, CorrelationState& state, bool first);

	int windowSize = 512;
	int hopLength = 4;
	int sampleRate = 48000.0;
	float pitch = 0.f;
	int startSample = 0;
	float correlationThreshold = 0.95;
	bool correlationStatus = false;
	bool hoppingEnabled = true;

	std::vector<float> rmsTransient;
};