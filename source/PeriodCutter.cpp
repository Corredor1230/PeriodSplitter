#include "PeriodCutter.h"

PeriodCutter::PeriodCutter(const Sitrano::AnalysisUnit& ana, Sitrano::Results& res, float sizeInMs, float hopInMs) :
	unit(ana), res(res), pitchDetector(res.pitch)
{
	sampleRate = unit.sampleRate;
}

void PeriodCutter::initialize(float newPitch, float threshold)
{
	setPitch(newPitch);
	int rmsLength = (int)(sampleRate / pitch) / 2;
	rmsTransient.resize(rmsLength);
	startSample = findStartTransient(0, rmsTransient, 1.f, 1.f, 3.0, threshold);
	windowSize = (int)(((float)sampleRate / pitch) + 1);
	allZeroes = findZeroCrossings(unit.soundFile, startSample);
	zeroList = findPeriodSamples(unit.soundFile, startSample, 50.0, pitch);
}

void PeriodCutter::setWindowSize(int sizeInSamples)
{
	windowSize = sizeInSamples;
}

void PeriodCutter::setWindowSize(float sizeInMs)
{
	windowSize = (sampleRate / 1000.0) * sizeInMs;
	for(int i = 0; i < 15; i++)
	{ 
		if (windowSize > pow(2, i) && windowSize < pow(2, i + 1))
			windowSize = pow(2, i + 1);
	}
}

void PeriodCutter::setThreshold(float thresh)
{
	correlationThreshold = thresh;
}

void PeriodCutter::setPitch(float newPitch)
{
	pitch = newPitch;
}

void PeriodCutter::enableHopping(bool hopping)
{
	hoppingEnabled = hopping;
}

void PeriodCutter::setHopLength(float hopInMs)
{
	hopLength = (int)(((float)sampleRate / 1000.0) * hopInMs);
}

void PeriodCutter::setHopLength(int hopInSamples)
{
	hopLength = hopInSamples;
}

void PeriodCutter::printZeroList(std::string& filename)
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

std::vector<int> PeriodCutter::getCorrelationPeaks()
{
	return peakList;
}

std::vector<int> PeriodCutter::getCorrelationZeroes()
{
	return zeroList;
}

int PeriodCutter::findStartTransient(int startSample, std::vector<float>& tVector, int rmsSize,
	int rmsHopLength, float transFactor, float transThreshold)
{
	std::vector<float> rmsList;
	float factor = transFactor;
	float threshold = transThreshold;
	bool transientFound = false;
	int tInitSample = 0;

	for (int samp = startSample; samp < unit.soundFile.size(); samp+=rmsHopLength)
	{
		float windowSum = 0.f;
		std::vector<float> tempRms(rmsSize);

		for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < unit.soundFile.size(); rmsSamp++)
		{
			float x = unit.soundFile[samp + rmsSamp];
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
	return findPreviousZero(unit.soundFile, peakSample);
}

int PeriodCutter::findStartTransient(int startSample, std::vector<float>& tVector, float rmsSizeInMs,
	float rmsHopLengthInRatio, float transFactor, float transThreshold)
{
	std::vector<float> rmsList;
	float factor = transFactor;
	float threshold = transThreshold;
	bool transientFound = false;
	int tInitSample = 0;
	int rmsSize = (rmsSizeInMs / 1000.f) * sampleRate;
	int rmsHopLength = (int)((float)rmsSize * rmsHopLengthInRatio);

	for (int samp = startSample; samp < unit.soundFile.size(); samp += rmsHopLength)
	{
		float windowSum = 0.f;
		std::vector<float> tempRms(rmsSize);


		for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < unit.soundFile.size(); rmsSamp++)
		{
			float x = unit.soundFile[samp + rmsSamp];
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
	return findPreviousZero(unit.soundFile, peakSample);
}

int PeriodCutter::findPeak(int startSample, std::vector<float> transientVector)
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

int PeriodCutter::findPreviousZero(const std::vector<float>& signal, int startSample)
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

int PeriodCutter::findNextZero(const std::vector<float>& signal, int startSample)
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

int PeriodCutter::findNearestZero(const std::vector<float>& signal, int startSample)
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

float PeriodCutter::findPeakValue(const std::vector<float>& signal, bool useAbsolute)
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

int PeriodCutter::findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute)
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

int PeriodCutter::findPeakSample(const std::vector<float>& signal, bool useAbsolute)
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

float PeriodCutter::findMode(const std::vector<float>& data, float threshold = 1.0f,
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

std::vector<int> PeriodCutter::findPeriodSamples(const std::vector<float>& signal, int startSample,
	float msOffset, float inPitch, float correlationThreshold)
{
	int expectedPeriodLength = static_cast<int>(sampleRate / inPitch);
	int expectedNumPeriods = static_cast<int>(signal.size() / expectedPeriodLength);
	int initialSample = startSample + static_cast<int>(sampleRate * (msOffset / 1000.0));
	int peakSample = findPeakSample(signal, initialSample, initialSample + windowSize, true);
	int periodStart = findNearestZero(signal, peakSample);
	std::vector<int> zeroCrossings = findZeroCrossings(signal, startSample);
	int margin = std::max<int>(4, std::ceil<int>((float)expectedPeriodLength * 0.02));

	std::vector<float> window(windowSize);
	for (int i = 0; i < windowSize; i++)
		window[i] = (periodStart + i < signal.size()) ? signal[periodStart + i] : 0.f;
	std::vector<float> backWindow = window;

	// Precompute window norm
	float squareA = 0.f;
	for (float w : window) squareA += w * w;
	if (squareA == 0.f) return {}; // Avoid divide by zero
	float normA = std::sqrt(squareA);

	std::deque<int> forwardList;
	std::deque<int> backwardList;
	forwardList.push_back(periodStart);
	backwardList.push_back(periodStart);

	auto forwardSearch = [&]() {
		std::vector<float> corrVals;
		CorrelationState state;
		int filePos = periodStart;
		bool first = true;
		int lastStart = filePos;

		while (filePos < signal.size()) {
			if (first) {
				lastStart = filePos;
				for (int i = 0; i < windowSize; i++)
					window[i] = (lastStart + i < signal.size()) ? signal[lastStart + i] : 0.f;
				for (float w : window) squareA += w * w;
				//if (squareA == 0.f) squareA = 0.0001; // Avoid divide by zero
				float normA = std::sqrt(squareA);
				first = false;
				filePos += (expectedPeriodLength - margin - 2);
			}

			float corr = signalCorrelationRolling(window, squareA, signal, filePos, state, false);
			if (corr > correlationThreshold) {
				corrVals.push_back(corr);
			}
			else if (filePos > lastStart + expectedPeriodLength + 3) {
				int nz = findNearestZeroCached(zeroCrossings, lastStart + expectedPeriodLength);
				if (nz == forwardList.back()) { filePos++; continue; }
				if (nz - forwardList.back() >= expectedPeriodLength - margin &&
					nz - forwardList.back() <= expectedPeriodLength + margin)
				{
					forwardList.push_back(nz);
					filePos = nz;
				}
				else
				{
					forwardList.push_back(lastStart + expectedPeriodLength);
					filePos = lastStart + expectedPeriodLength;
				}
				first = true;
			}
			else {
				filePos++;
				continue;
			}

			if (corrVals.size() < 3) continue;
			int lc = corrVals.size() - 1;
			bool isPeak = (corrVals[lc - 1] > corrVals[lc]) && (corrVals[lc - 1] > corrVals[lc - 2]);
			bool inRange = filePos > forwardList.back() + expectedPeriodLength - margin &&
				filePos < forwardList.back() + expectedPeriodLength + margin;

			if (isPeak /*&& inRange*/) {
				int nz = findNearestZeroCached(zeroCrossings, filePos);
				if (nz == forwardList.back()) { filePos++; continue; }
				int currentWindowSize = nz - forwardList.back();
				if (currentWindowSize >= (expectedPeriodLength - margin)/* &&
					currentWindowSize <= (expectedPeriodLength + margin)*/)
				{
					forwardList.push_back(nz);
					filePos = nz;
				}
				else
				{
					int temporaryNewStart = forwardList.back() + expectedPeriodLength;
					forwardList.push_back(temporaryNewStart);
					filePos = temporaryNewStart;
				}
				first = true;
			}
			if (forwardList.size() > expectedNumPeriods * 1.1) break;
			filePos++;
		}
	};

	auto backwardSearch = [&]() {
		std::vector<float> corrVals;
		CorrelationState state;
		int filePos = periodStart;
		bool first = true;
		int lastStart = filePos;

		while (filePos > (startSample - windowSize / 2)) {
			if (first) {
				first = false;
				lastStart = filePos;
				for (int i = 0; i < windowSize; i++)
					backWindow[i] = (lastStart + i) > 1 ? signal[lastStart + i] : 0.f;
				for (float w : backWindow) squareA += w * w;
				//if (squareA == 0.f) squareA = 0.0001; // Avoid divide by zero
				float normA = std::sqrt(squareA);
				filePos -= (expectedPeriodLength - margin);
				corrVals.push_back(0.0);
			}

			//float corr = signalCorrelationRolling(backWindow, squareA, signal, filePos, state, false);
			float corr = signalCorrelation(backWindow, signal, filePos);
			corrVals.push_back(corr);

			if (corrVals.size() < 3)
			{
				filePos--;
				continue;
			}
			int lc = corrVals.size() - 1;

			bool aboveThreshold = (
				(corrVals[lc] > correlationThreshold) ||
				(corrVals[lc - 1] > correlationThreshold) ||
				(corrVals[lc - 2] > correlationThreshold)
				);

			bool isPeak = (
				(corrVals[lc - 1] > corrVals[lc]) && 
				(corrVals[lc - 1] > corrVals[lc - 2])
				);

			bool inRange = (
				(filePos < backwardList.front() - expectedPeriodLength + margin) &&
				(filePos > backwardList.front() - expectedPeriodLength - margin)
				);

			if (isPeak  /*&& inRange*/ && aboveThreshold) {
				int nz = findNearestZeroCached(zeroCrossings, filePos);
				int currentWindowSize = backwardList.front() - nz;
				if (currentWindowSize <= expectedPeriodLength + margin &&
					currentWindowSize >= expectedPeriodLength - margin)
				{
					backwardList.push_front(nz);
					filePos = nz + 1;
				}
				else
				{
					backwardList.push_front(lastStart - expectedPeriodLength);
					filePos = lastStart - expectedPeriodLength;
				}
				first = true;
			}
			filePos--;
		}
	};

	// Run forward and backward search in parallel
	auto forwardFuture = std::async(std::launch::async, forwardSearch);
	auto backwardFuture = std::async(std::launch::async, backwardSearch);
	forwardFuture.get();
	backwardFuture.get();

	// Merge both lists
	std::vector<int> periodZeroes;
	periodZeroes.insert(periodZeroes.end(), backwardList.begin(), backwardList.end());
	periodZeroes.insert(periodZeroes.end(), ++forwardList.begin(), forwardList.end());

	return periodZeroes;
}


std::vector<int> PeriodCutter::findZeroCrossings(const std::vector<float>& signal, int initSample)
{
	std::vector<int> zeroCrossings;
	for (int i = initSample; i < signal.size(); ++i)
	{
		if ((signal[i - 1] < 0 && signal[i] >= 0) ||
			(signal[i - 1] > 0 && signal[i] <= 0))
		{
			zeroCrossings.push_back(i);
		}
	}

	return zeroCrossings;
}

int PeriodCutter::findNearestZeroCached(const std::vector<int>& zeroCrossings, int sample)
{
	auto it = std::lower_bound(zeroCrossings.begin(), zeroCrossings.end(), sample);
	if (it == zeroCrossings.end()) return zeroCrossings.back();
	if (it == zeroCrossings.begin()) return *it;
	int after = *it;
	int before = *(it - 1);
	return (std::abs(after - sample) < std::abs(before - sample)) ? after : before;
}

float PeriodCutter::findPitch(const std::vector<float>& signal)
{
	std::vector<float> pitches = pitchDetector.feed(signal);

	pitch = findMode(pitches);

	std::cout << "Current pitch: " << pitch << "Hz\n";

	return pitch;
}

float PeriodCutter::signalCorrelation(std::vector<float>& window, const std::vector<float>& signal, int startSample)
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

float PeriodCutter::signalCorrelationRolling(const std::vector<float>& window, float squareA,
	const std::vector<float>& signal, int startSample, CorrelationState& state, bool first)
{
	state.numerator = 0.f;
	state.squareB = 0.f;

	for (int i = 0; i < window.size(); ++i)
	{
		float a = window[i];
		float b = (i + startSample < signal.size()) ? signal[i + startSample] : 0.f;
		state.numerator += a * b;
		state.squareB += b * b;
	}

	float denominator = std::sqrt(squareA * state.squareB);
	return (denominator != 0.f) ? state.numerator / denominator : 0.f;
}