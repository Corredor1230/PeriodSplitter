#include "Correlator.h"

Correlator::Correlator(std::vector<float>& file, SF_INFO& info, float sizeInMs, float hopInMs) :
	audioFile(file), sfInfo(info), pitchDetector(sfInfo.samplerate)
{
	sampleRate = sfInfo.samplerate;

	if (sfInfo.channels > 1)
		separateChannels(audioFile);
	correlation.resize(audioFile.size() + windowSize);

	findPitch();

	windowSize = (int)(((float)sampleRate / pitch) + 0.1 * (float)sampleRate / (pitch));

	int rmsLength = (int)(sampleRate / pitch) / 2;
	rmsTransient.resize(rmsLength);

	startSample = findStartTransient(0, rmsTransient, 1.f, 1.f, 3.0, 0.1);
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

void Correlator::moveToNextWindow(std::vector<float>& window, int currentSamp)
{
	for (int samp = 0; samp < windowSize; samp++)
	{
		if (samp < window.size())
		{
			if (currentSamp + windowSize + samp < audioFile.size())
				window[samp] = audioFile[currentSamp + samp];
			else
				window[samp] = 0.f;
		}
	}
}

void Correlator::printCorrelation(std::string& filename)
{
	if (filename.empty())
	{
		std::ofstream outfile("AudioCorrelation.txt");
		for (int i = 0; i < correlation.size(); i++)
		{
			if(correlation[i] != 0.0)
				outfile << "Sample: " << i << " | Correlation: " << correlation[i] << std::endl;
		}
	}
	else
	{
		std::ofstream outfile(filename + ".txt");
		for (int i = 0; i < correlation.size(); i++)
		{
			if(correlation[i] != 0.0)
				outfile << "Sample: " << i << " | Correlation: " << correlation[i] << std::endl;
		}
	}
}

void Correlator::printCorrelationPeak(std::string& filename)
{
	if (filename.empty())
	{
		std::ofstream outfile("AudioPeakCorrelation.txt");
		for (int i = 0; i < correlation.size(); i++)
		{
			outfile << "Sample: " << i << " | Correlation: " << correlation[i] << std::endl;
		}
	}
	else
	{
		std::ofstream outfile(filename + ".txt");
		for (int i = 0; i < peakList.size(); i++)
		{
			outfile << "Sample: " << peakList[i] << " | Correlation: " << correlation[peakList[i]] << std::endl;
		}
	}
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

void Correlator::calculateCorrelation(int startSample)
{
	int numWindows = (int)std::ceil((float)sfInfo.frames / (float)windowSize);
	std::vector<float> corrWindow(windowSize);
	int cumulativeGap = 0;
	bool firstPeakFound = false;
	int lastSample = 0;
	int offset = 0;
	for (int windowNum = 0; windowNum < numWindows; windowNum++)
	{
		float peakValue = 0.f;
		bool correlationInWindow = false;
		int numberOfCorrelations = 0;

		//Load initial window
		if(windowNum == 0)
			moveToNextWindow(corrWindow, startSample + (windowNum * windowSize));

		//Set initial sample from original signal
		for (int initSamp = 0; initSamp < windowSize; initSamp++)
		{
			float numerator	= 0.f;
			float denominator = 0.f;
			float squareA = 0.f;
			float squareB = 0.f;
			int fileSample = 0;

			if (zeroList.empty())
				fileSample = initSamp - cumulativeGap + startSample + (windowNum * windowSize);
			else
				fileSample = zeroList.back() + initSamp + startSample;
			lastSample = fileSample;

			//Iterate through all window samples 
			//and calculate correlation
			for (int samp = 0; samp < windowSize; samp++)
			{
				float a = 0.f;

				//Check we're within range
				if (fileSample + samp < audioFile.size())
					a = audioFile[fileSample + samp];
				else
					a = 0.f;

				float b = corrWindow[samp];
				numerator += a * b;
				squareA += a * a;
				squareB += b * b;
			}
			denominator = std::sqrt(squareA * squareB);
			float corr = 0.f;

			if (denominator != 0)
				corr = numerator / denominator;
			else
				corr = 0.f;

			if (fileSample < correlation.size())
				correlation[fileSample] = corr;
			else
				correlation.push_back(corr);

			//Skips to next loop for initial loops
			if (initSamp < hopLength)
				continue;

			//Skips several loops if correlation is far from threshold
			//always subtract 1 because initSamp will increase 1 regardless
			if (initSamp < (audioFile.size() - hopLength * 2) && corr < 0.f)
			{
				initSamp += hopLength * 2 - 1;
				continue;
			}
			else if (initSamp < (audioFile.size() - hopLength) && corr < correlationThreshold)
			{
				initSamp += hopLength - 1;
				continue;
			}

			bool isPeak = correlation[fileSample - 1] > correlation[fileSample - 2]
				&& correlation[fileSample - 1] > correlation[fileSample];

			//Add peaks to peaklist
			if (isPeak && corr>correlationThreshold)
			{
				correlationInWindow = true;
				numberOfCorrelations++;
				peakList.push_back(fileSample - 1);
			}
		}

		//Checks that we've already found the first window's peak
		if (correlationInWindow && firstPeakFound)
		{
			int lastZero = 0;
			for (int i = 0; i < numberOfCorrelations; i++)
			{
				zeroList.push_back(0.0);
			}
			for (int i = 0; i < numberOfCorrelations; i++)
			{
				int zeroSample = findPreviousZero(peakList[peakList.size() - i - 1]);
				if (zeroSample != 0)
				{
					zeroList[zeroList.size() - i - 1] = zeroSample;
					if (i == 0)
						lastZero = zeroList.back();
				}
				else
					moveToNextWindow(corrWindow, startSample + ((windowNum + 1) * windowSize));
			}
			moveToNextWindow(corrWindow, lastZero);
			cumulativeGap += lastSample - lastZero;
		}

		else if (correlationInWindow && !firstPeakFound)
		{
			int zeroSample = 0;
			int peakSample = findPeakSample(corrWindow, true);
			zeroSample = findPreviousZero(peakSample + windowNum * windowSize);
			moveToNextWindow(corrWindow, zeroSample);
			cumulativeGap += lastSample - zeroSample;
			firstPeakFound = true;
		}

		else
		{
			moveToNextWindow(corrWindow, startSample + ((windowNum + 1) * windowSize));
		}
		std::cout << "Processed: " << (float)windowNum * 100 / (float)numWindows << "%\n";

		if (cumulativeGap > windowSize)
		{
			windowNum--;
			cumulativeGap = cumulativeGap - windowSize;
		}
	}

	//errorCorrection(audioFile, peakList);

	correlationStatus = true;
}

//void Correlator::calculateCorrelation(int startSample)
//{
//	int transientSample = findStartTransient(0, 128, 128);
//	for (int oSamp = transientSample; oSamp < audioFile.size(); oSamp++)
//	{
//		float numerator = 0.f;
//		float denominator = 0.f;
//		float currentCorrelation = 0.f;
//		bool validCorrelation = false;
//		float sqA = 0.f;
//		float sqB = 0.f;
//
//		for (int wSamp = 0; wSamp < windowSize; wSamp++)
//		{
//			float a = 0.f;
//			float b = 0.f;
//
//			if (wSamp + oSamp < audioFile.size())
//				a = audioFile[oSamp + wSamp];
//			else
//				a = 0.f;
//
//			b = audioFile[oSamp];
//
//			numerator += a * b;
//			sqA += a * a;
//			sqB += b * b;
//		}
//
//		denominator = std::sqrt(sqA * sqB);
//
//		//Calculate and log correlation
//		currentCorrelation = numerator / denominator;
//		correlation[oSamp] = currentCorrelation;
//
//		//Check correlation value and skip if there's no valid correlation
//		validCorrelation = currentCorrelation > correlationThreshold ? true : false;
//		if (!validCorrelation || oSamp < 2)
//		{
//			//Jump by hop length if correlation is larger than 0 but smaller than thresh
//			//jump by double if it's below 0
//			currentCorrelation > 0.f ? oSamp += hopLength / 2: oSamp += hopLength;
//			continue;
//		}
//
//		bool isPeak =
//			correlation[oSamp - 1] > correlation[oSamp - 2] &&
//			correlation[oSamp - 1] > correlation[oSamp];
//
//		if (validCorrelation && isPeak)
//		{
//			peakList.push_back(oSamp - 1);
//		}
//
//
//
//
//	}
//}

std::vector<int> Correlator::getCorrelationPeaks()
{
	return peakList;
}

std::vector<float> Correlator::getCorrelationVector()
{
	return correlation;
}

std::vector<int> Correlator::getCorrelationZeroes()
{
	return zeroList;
}

bool Correlator::hasBeenCalculated()
{
	return correlationStatus;
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
	return findPreviousZero(peakSample);
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
	return findPreviousZero(peakSample);
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

void Correlator::separateChannels(std::vector<float>& audioFile)
{
	int numChans = sfInfo.channels;
	int numSamps = audioFile.size();

	for (int chan = 0; chan < numChans; chan++)
	{
		std::vector<float> chanVector(numSamps / numChans);

		for (int samp = 0; samp < chanVector.size(); samp++)
		{
			chanVector[samp] = audioFile[(samp * numChans) + chan];
		}
		multiChannel.push_back(chanVector);
	}
}

int Correlator::findPreviousZero(int startSample)
{
	int zeroSample = 0;
	for (int i = startSample; i > 2; i--)
	{
		bool zeroFound = false;

		zeroFound = (audioFile[i] > 0.0 && audioFile[i - 2] < 0.0)
			  || (audioFile[i] < 0.0 && audioFile[i - 2] > 0.0);

		if (zeroFound)
		{
			if (std::abs(audioFile[i - 2]) < std::abs(audioFile[i - 1]))
				zeroSample = i - 2;
			else
				zeroSample = i - 1;
			break;
		}
	}

	return zeroSample;
}

int Correlator::findNearestZero(int startSample)
{
	int prevZero = 0;
	int distanceA = 0;
	int nextZero = 0;
	int distanceB = 0;

	for (int i = startSample; i > 2 && i < audioFile.size(); i--)
	{
		bool zeroFound = false;

		zeroFound = (audioFile[i] > 0.0 && audioFile[i - 2] < 0.0)
			|| (audioFile[i] < 0.0 && audioFile[i - 2] > 0.0);

		if (zeroFound)
		{
			if (std::abs(audioFile[i - 1]) < std::abs(audioFile[i - 2]))
				prevZero = i - 1;
			else
				prevZero = i - 2;
			distanceA = std::abs(startSample - prevZero);
			break;
		}
	}

	for (int i = startSample; i < audioFile.size() - 2; i++)
	{
		bool zeroFound = false;

		zeroFound = (audioFile[i] > 0.0 && audioFile[i + 2] < 0.0)
			|| (audioFile[i] < 0.0 && audioFile[i + 2] > 0.0);

		if (zeroFound)
		{
			if (std::abs(audioFile[i + 1]) < std::abs(audioFile[i + 2]))
				prevZero = i + 1;
			else
				prevZero = i + 2;
			nextZero = i + 1;
			distanceB = std::abs(startSample - nextZero);
			break;
		}
	}

	return std::min(distanceA, distanceB);
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

void Correlator::findPitch()
{
	std::vector<float> pitches = pitchDetector.feed(audioFile);

	pitch = findMode(pitches);

	std::cout << "Current pitch: " << pitch << "Hz\n";
}

void Correlator::errorCorrection(std::vector<float>& signal, std::vector<int>& peaks, float errorThreshold)
{
	int totalNumber = 0;
	for (int i = 0; i < peaks.size() - 1; i++)
	{
		totalNumber += peaks[i + 1] - peaks[i];
	}
	int averageSize = totalNumber / peaks.size();
	int errorMargin = (int)(averageSize - ((float)averageSize * errorThreshold));

	std::vector<int> errorLoops;
	for (int i = 0; i < peaks.size() - 1; i++)
	{
		int loopSize = peaks[i + 1] - peaks[i];
		if (loopSize < averageSize - errorMargin ||
			loopSize > averageSize + errorMargin)
		{
			errorLoops.push_back(i);
			std::cout << "Loop No. " << i << " || Size: " << loopSize << '\n';
		}
	}

	std::cout << "Total number of errors: " << errorLoops.size() << '\n';

	for (int i = 0; i < errorLoops.size(); i++)
	{

	}


}