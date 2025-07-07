#include "Correlator.h"

Correlator::Correlator(std::vector<float>& file, SF_INFO& info, float sizeInMs, float hopInMs) :
	audioFile(file), sfInfo(info), pitchDetector(sfInfo.samplerate)
{
	sampleRate = sfInfo.samplerate;

	pitch = findPitch(audioFile);
	int rmsLength = (int)(sampleRate / pitch) / 2;
	rmsTransient.resize(rmsLength);

	startSample = findStartTransient(0, rmsTransient, 1.f, 1.f, 3.0, 0.1);

	windowSize = (int)(((float)sampleRate / pitch) + 0.05 * (float)sampleRate / (pitch));
	zeroList = findPeriodSamples(file, startSample, 50.0, pitch);
}

void Correlator::initialize(int sr, int sizeInSamples)
{
	sampleRate = sr;
	setWindowSize(sizeInSamples);
}

void Correlator::initialize(int sr, float sizeInMs)
{
	sampleRate = sr;
	setWindowSize(sizeInMs);
}

void Correlator::setWindowSize(int sizeInSamples)
{
	windowSize = sizeInSamples;
}

void Correlator::setWindowSize(float sizeInMs)
{
	windowSize = (sampleRate / 1000.0) * sizeInMs;
	for(int i = 0; i < 15; i++)
	{ 
		if (windowSize > pow(2, i) && windowSize < pow(2, i + 1))
			windowSize = pow(2, i + 1);
	}
}

void Correlator::setThreshold(float thresh)
{
	correlationThreshold = thresh;
}

void Correlator::enableHopping(bool hopping)
{
	hoppingEnabled = hopping;
}

void Correlator::setHopLength(float hopInMs)
{
	hopLength = (int)(((float)sampleRate / 1000.0) * hopInMs);
}

void Correlator::setHopLength(int hopInSamples)
{
	hopLength = hopInSamples;
}

void Correlator::printZeroList(std::string& filename)
{
	if (filename.empty())
	{
		std::ofstream outfile("AudioZeroList.txt");
		for (int i = 0; i < zeroList.size(); i++)
		{
			outfile << "Sample: " << zeroList[i] << std::endl;
		}
	}
	else
	{
		std::ofstream outfile(filename + ".txt");
		for (int i = 0; i < zeroList.size(); i++)
		{
			outfile << "Sample: " << zeroList[i] << std::endl;
		}
	}
}

std::vector<int> Correlator::getCorrelationPeaks()
{
	return peakList;
}

std::vector<int> Correlator::getCorrelationZeroes()
{
	return zeroList;
}

int Correlator::findStartTransient(int startSample, std::vector<float>& tVector, int rmsSize,
	int rmsHopLength, float transFactor, float transThreshold)
{
	std::vector<float> rmsList;
	float factor = transFactor;
	float threshold = transThreshold;
	bool transientFound = false;
	int tInitSample = 0;

	for (int samp = startSample; samp < audioFile.size(); samp+=rmsHopLength)
	{
		float windowSum = 0.f;
		std::vector<float> tempRms(rmsSize);

		for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < audioFile.size(); rmsSamp++)
		{
			float x = audioFile[samp + rmsSamp];
			windowSum += (x * x);
			tempRms[rmsSamp] = x;
		}

		float rms = std::sqrt(windowSum / rmsSize);
		if (rmsList.size() < 2)
		{
			rmsList.push_back(rms);
			continue;
		}

		float rmsRatio = rmsList.back() == 0.f ? 1.f : rms / rmsList.back();
		rmsList.push_back(rms);

		if (rmsRatio > factor && rms > threshold)
		{
			transientFound = true;
			if (rmsTransient.size() != tempRms.size())
				rmsTransient.resize(tempRms.size());
			rmsTransient = tempRms;
			tInitSample = samp;
			break;
		}
	}
	int rmsPeakSample = findPeak(0, rmsTransient);
	int peakSample = tInitSample + rmsPeakSample;
	return findPreviousZero(audioFile, peakSample);
}

int Correlator::findStartTransient(int startSample, std::vector<float>& tVector, float rmsSizeInMs,
	float rmsHopLengthInRatio, float transFactor, float transThreshold)
{
	std::vector<float> rmsList;
	float factor = transFactor;
	float threshold = transThreshold;
	bool transientFound = false;
	int tInitSample = 0;
	int rmsSize = (rmsSizeInMs / 1000.f) * sampleRate;
	int rmsHopLength = (int)((float)rmsSize * rmsHopLengthInRatio);

	for (int samp = startSample; samp < audioFile.size(); samp += rmsHopLength)
	{
		float windowSum = 0.f;
		std::vector<float> tempRms(rmsSize);


		for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < audioFile.size(); rmsSamp++)
		{
			float x = audioFile[samp + rmsSamp];
			windowSum += (x * x);
			tempRms[rmsSamp] = x;
		}

		float rms = std::sqrt(windowSum / rmsSize);
		if (rmsList.size() < 2)
		{
			rmsList.push_back(rms);
			continue;
		}

		float rmsRatio = rmsList.back() == 0.f ? 1.f : rms / rmsList.back();
		rmsList.push_back(rms);

		if (rmsRatio > factor && rms > threshold)
		{
			transientFound = true;
			if (rmsTransient.size() != tempRms.size())
				rmsTransient.resize(tempRms.size());
			rmsTransient = tempRms;
			tInitSample = samp;
			break;
		}
	}
	int rmsPeakSample = findPeak(0, rmsTransient);
	int peakSample = tInitSample + rmsPeakSample;
	return findPreviousZero(audioFile, peakSample);
}

int Correlator::findPeak(int startSample, std::vector<float> transientVector)
{
	std::vector<float> absTransient(transientVector.size());

	for (int i = 0; i < absTransient.size(); i++)
	{
		absTransient[i] = std::abs(transientVector[i]);
	}

	int peakSample = std::distance(absTransient.begin(),
		std::max_element(absTransient.begin(), absTransient.end()));
	return peakSample;
}

int Correlator::findPreviousZero(std::vector<float>& signal, int startSample)
{
	int zeroSample = 0;
	for (int i = startSample; i > 2; i--)
	{
		bool zeroFound = false;

		zeroFound = (signal[i] >= 0.0 && signal[i - 2] <= 0.0)
			  || (signal[i] <= 0.0 && signal[i - 2] >= 0.0);

		if (zeroFound)
		{
			if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
				zeroSample = i - 2;
			else
				zeroSample = i - 1;
			break;
		}
	}

	return zeroSample;
}

int Correlator::findNextZero(std::vector<float>& signal, int startSample)
{
	int zeroSample = 0;
	for (int i = startSample; i < signal.size(); i++)
	{
		bool zeroFound = false;

		zeroFound = (signal[i] >= 0.0 && signal[i - 2] < 0.0)
			|| (signal[i] <= 0.0 && signal[i - 2] > 0.0);

		if (zeroFound)
		{
			if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
				zeroSample = i - 2;
			else
				zeroSample = i - 1;
			break;
		}
	}

	return zeroSample;
}

int Correlator::findNearestZero(std::vector<float>& signal, int startSample)
{
	int prevZero = findPreviousZero(signal, startSample);
	int nextZero = findNextZero(signal, startSample);

	int distancePrev = std::abs(startSample - prevZero);
	int distanceNext = std::abs(startSample - nextZero);

	if (distancePrev < distanceNext)
		return prevZero;
	else if (distanceNext < distancePrev)
		return nextZero;
	else
		return prevZero;
}

float Correlator::findPeakValue(std::vector<float>& signal, bool useAbsolute)
{
	float peakValue = 0.f;
	for (int i = 0; i < signal.size(); i++)
	{
		if (useAbsolute)
		{
			if (std::abs(signal[i]) > peakValue)
				peakValue = std::abs(signal[i]);
		}
		else
		{
			if (signal[i] > peakValue)
				peakValue = signal[i];
		}
	}
	return peakValue;
}

int Correlator::findPeakSample(std::vector<float>& signal, int startSample, int endSample, bool useAbsolute)
{
	int peakSample = 0;
	float peakValue = 0.f;
	int size = endSample - startSample;
	for (int i = 0; i < size; i++)
	{
		int sample = startSample + i;
		if (useAbsolute)
		{
			if (std::abs(signal[sample]) > peakValue)
			{
				peakValue = std::abs(signal[sample]);
				peakSample = sample;
			}
		}
		else
		{
			if (signal[sample] > peakValue)
			{
				peakValue = signal[sample];
				peakSample = sample;
			}
		}
	}

	return peakSample;
}

int Correlator::findPeakSample(std::vector<float>& signal, bool useAbsolute)
{
	int peakSample = 0;
	float peakValue = 0.f;
	for (int i = 0; i < signal.size(); i++)
	{
		if (useAbsolute)
		{
			if (std::abs(signal[i]) > peakValue)
			{
				peakValue = std::abs(signal[i]);
				peakSample = i;
			}
		}
		else
		{
			if (signal[i] > peakValue)
			{
				peakValue = signal[i];
				peakSample = i;
			}
		}
	}
	return peakSample;
}

float Correlator::findMode(const std::vector<float>& data, float threshold = 1.0f,
	float minFreq = 60.0, float maxFreq = 1200.0) {
	if (data.empty()) return NAN;

	// Step 1: Sort the data
	std::vector<float> sortedData = data;
	std::sort(sortedData.begin(), sortedData.end());

	// Step 2: Slide a window and track the group with the highest count
	int maxCount = 0;
	float modeSum = 0.0f;
	int modeCount = 0;

	for (size_t i = 0; i < sortedData.size(); ++i) {
		float ref = sortedData[i];
		float sum = 0.0f;
		int count = 0;

		if (ref < minFreq || ref > maxFreq)
			continue;

		// Include all values within [ref, ref + threshold)
		for (size_t j = i; j < sortedData.size() && sortedData[j] - ref < threshold; ++j) {
			sum += sortedData[j];
			++count;
		}

		if (count > maxCount) {
			maxCount = count;
			modeSum = sum;
			modeCount = count;
		}
	}

	float mode = modeSum / modeCount;

	return modeCount > 0 ? mode : NAN;
}

std::vector<int> Correlator::findPeriodSamples(std::vector<float>& signal, int startSample,
	float msOffset, float inPitch, float correlationThreshold)
{
	int expectedPeriodLength = (int)(sampleRate / inPitch);
	int expectedNumPeriods = (int)(signal.size() / expectedPeriodLength);
	int initialSample = startSample + (int)(sampleRate * (msOffset / 1000.0));
	int peakSample = findPeakSample(signal, initialSample, initialSample + windowSize, true);
	int periodStart = findPreviousZero(signal, peakSample);
	float threshold = correlationThreshold;
	std::deque<int> periodList;
	periodList.push_back(periodStart);

	std::vector<float> window(windowSize);
	std::vector<float> correlationValues;

	int hopSize = (int)(sampleRate / pitch);
	bool firstInPeriod = true;
	int lastStart = 0;

	
	/* This loop will go over all the samples in the full audio file.
	* The inner loop for the window is within signalCorrelation() */
	for (int filePos = periodStart; filePos < signal.size(); filePos++)
	{
		if (firstInPeriod)
		{
			lastStart = filePos;
			//Refreshes the window
			for (int i = 0; i < windowSize; i++)
			{
				if (filePos + i >= signal.size())
					window[i] = 0.0;
				else
					window[i] = signal[filePos + i];
			}

			firstInPeriod = false;

			//At very high frequencies this will just be 0
			//But at low frequencies it saves a lot of time.
			filePos += (expectedPeriodLength * 3) / 4;
		}

		float currentCorrelation = signalCorrelation(window, signal, filePos);

		//Verifies that we have a significant correlation
		if (currentCorrelation > correlationThreshold)
		{
			correlationValues.push_back(currentCorrelation);
		}
		else if (filePos > (int)(lastStart + ((float)expectedPeriodLength * 1.1)))
		{
			int nearestZero = findNearestZero(signal, lastStart + expectedPeriodLength + 2);
			if (nearestZero == periodList.back())
				continue;
			periodList.push_back(nearestZero);
			filePos = nearestZero;
			firstInPeriod = true;
		}
		else
			continue;

		//We cannot check for any correlation peaks with less than 3 values
		//So we just save them and skip to the next loop.
		if (correlationValues.size() < 3)
			continue;

		//This checks that the correlation value is a peak
		int lastCorr = correlationValues.size() - 1;
		bool isCorrelationPeak = (correlationValues[lastCorr - 1] > correlationValues[lastCorr])
			&& (correlationValues[lastCorr - 1] > correlationValues[lastCorr - 2]);
		//This checks that we're within 10% of the expected period size
		bool isWithinPeriodRange = filePos > (int)(((float)periodList.back() + (expectedPeriodLength * 0.9))) 
			&& filePos < (int)(((float)periodList.back() + (expectedPeriodLength * 1.1)));

		if (isCorrelationPeak && isWithinPeriodRange)
		{
			int nearestZero = findNearestZero(signal, filePos);
			if (nearestZero == periodList.back())
				continue;
			periodList.push_back(nearestZero);
			filePos = nearestZero;
			firstInPeriod = true;
		}
		else
			continue;

		if (periodList.size() > (int)((float)expectedNumPeriods * 1.1))
			break;
	}

	firstInPeriod = true;
	correlationValues.clear();

	for (int filePos = periodStart; filePos > startSample; filePos--)
	{
		if (firstInPeriod)
		{
			//Refreshes the window
			for (int i = 0; i < windowSize; i++)
			{
				window[i] = signal[filePos + i];
			}

			firstInPeriod = false;

			//At very high frequencies this will just be 0
			//But at low frequencies it saves a lot of time.
			filePos -= (expectedPeriodLength * 3) / 4;
		}

		float currentCorrelation = signalCorrelation(window, signal, filePos);

		//Verifies that we have a significant correlation
		if (currentCorrelation > correlationThreshold)
		{
			correlationValues.push_back(currentCorrelation);
		}
		else
			continue;

		//We cannot check for any correlation peaks with less than 3 values
		//So we just save them and skip to the next loop.
		if (correlationValues.size() < 3)
			continue;

		//This checks that the correlation value is a peak
		int lastCorr = correlationValues.size() - 1;
		bool isCorrelationPeak = (correlationValues[lastCorr - 1] > correlationValues[lastCorr])
			&& (correlationValues[lastCorr - 1] > correlationValues[lastCorr - 2]);
		//This checks that we're within 10% of the expected period size
		bool isWithinPeriodRange = filePos < (int)(((float)periodList.front() - (expectedPeriodLength * 1.1)))
			&& filePos > (int)(((float)periodList.front() - (expectedPeriodLength * 0.9)));

		if (isCorrelationPeak && isWithinPeriodRange)
		{
			int nearestZero = findNearestZero(signal, filePos);
			periodList.push_front(nearestZero);
			filePos = nearestZero;
			firstInPeriod = true;
		}
		else
			continue;
	}

	std::vector<int> periodZeroes(periodList.size());

	for (int i = 0; i < periodList.size(); i++)
	{
		periodZeroes[i] = periodList[i];
	}

	return periodZeroes;
}

float Correlator::findPitch(const std::vector<float>& signal)
{
	std::vector<float> pitches = pitchDetector.feed(signal);

	pitch = findMode(pitches);

	std::cout << "Current pitch: " << pitch << "Hz\n";

	return pitch;
}

float Correlator::signalCorrelation(std::vector<float>& window, std::vector<float>& signal, int startSample)
{
	float squareA = 0.f;
	float squareB = 0.f;
	float denominator = 0.f;
	float numerator = 0.f;
	float correlationValue = 0.f;
	for (int i = 0; i < window.size(); i++)
	{
		float a = window[i];
		float b = 0.f;
		if ((i + startSample) < signal.size())
			b = signal[i + startSample];
		else
			b = 0.f;

		numerator += a * b;
		squareA += a * a;
		squareB += b * b;
	}
	denominator = std::sqrt(squareA * squareB);
	if (denominator != 0.f)
		correlationValue = numerator / denominator;
	else
		correlationValue = 0.f;

	return correlationValue;
}