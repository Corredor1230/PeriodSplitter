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

class Correlator
{
public:
	Correlator(std::vector<float>& file, SF_INFO& info, float sizeInMs = 50.f, float hopInMs = 0.25);
	~Correlator() {};

	void initialize(int sampleRate, int sizeInSamples);
	void initialize(int sampleRate, float sizeInMs);
	void setWindowSize(int sizeInSamples = 512);
	void setWindowSize(float sizeInMs);
	void setThreshold(float thresh);
	void enableHopping(bool hopping);
	void setHopLength(float hopInMs);
	void setHopLength(int hopInSamples);
	void printCorrelation(std::string& filename);
	void printCorrelationPeak(std::string& filename);
	void printZeroList(std::string& filename);
	void calculateCorrelation(int startSample);
	std::vector<int> getCorrelationPeaks();
	std::vector<float> getCorrelationVector();
	std::vector<int> getCorrelationZeroes();
	bool hasBeenCalculated();

	std::vector<std::vector<int>>	chopPositions;
	std::vector<std::vector<float>>	multiChannel;

private:

	std::vector<float>&	audioFile;
	std::vector<float>	correlation;
	std::vector<int>	peakList;
	std::vector<int>	zeroList;
	SF_INFO&			sfInfo;
	PyinCpp				pitchDetector;

	int findStartTransient(int startSample, std::vector<float>& vector, int rmsSize, 
		int rmsHopLength,float transFactor = 2.0, float transThresh = 0.2);
	int findStartTransient(int startSample, std::vector<float>& vector, float rmsSizeInMs,
		float rmsHopLengthInRatio, float transFactor = 2.0, float transThresh = 0.2);
	int findPeak(int startSample, std::vector<float> vector);

	void moveToNextWindow(std::vector<float>& window, int currentSamp);
	void separateChannels(std::vector<float>& file);
	int findPreviousZero(int startSample);
	int findNearestZero(int startSample);
	float findPeakValue(std::vector<float>& signal, bool useAbsolute = true);
	int findPeakSample(std::vector<float>& signal, int startSample, int endSample, bool useAbsolute = true);
	int findPeakSample(std::vector<float>& signal, bool useAbsolute = true);
	std::deque<int> findPeriodSamples(std::vector<float>& signal, int startSample, 
		float msOffset, float inPitch, float correlationThreshold = 0.8);
	float findMode(const std::vector<float>& vctr, float threshold, float minFreq, float maxFreq);
	void findPitch();
	void errorCorrection(std::vector<float>& signal, std::vector<int>& peaks, float errorThreshold = 0.95);
	float signalCorrelation(std::vector<float>& window, std::vector<float>& signal, int startSample);


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