#include"TransientRecognition.h"

Transient::Transient(
	const Sitrano::AnalysisUnit& unit,
    const Sitrano::TransientSettings& conf
) :
    tSettings(conf),
	sampleRate(unit.sampleRate),
	aud(unit.soundFile),
    factor(conf.transientFactor),
    threshold(conf.transientThreshold)
{
    if (!conf.useMs)
    {
        rmsSize = conf.rmsSampleSize;
        rmsHopLength = conf.rmsSampleHopLength;
    }
    else
    {
        rmsSize = (tSettings.transientRmsSizeMs / 1000.f) * sampleRate;
        rmsHopLength = std::max(1, (int)((float)rmsSize * tSettings.transientRmsHopRatio));
    }
}

Sitrano::SampleRange Transient::findStartTransient()
{
    Sitrano::SampleRange range;
    range.initSample = tSettings.tStartSample;

    std::vector<float> rmsList;
    bool transientFound = false;
    int tInitSample = 0;

    std::vector<float> rmsWindow(rmsSize); // This will hold the transient

    for (int samp = tSettings.tStartSample; samp < aud.size(); samp += rmsHopLength)
    {
        float windowSum = 0.f;
        std::vector<float> tempRms(rmsSize);

        for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < aud.size(); rmsSamp++)
        {
            float x = aud[samp + rmsSamp];
            windowSum += (x * x);
            tempRms[rmsSamp] = x;
        }

        float rms = std::sqrt(windowSum / rmsSize);
        if (rmsList.size() < 2) {
            rmsList.push_back(rms);
            continue;
        }

        float rmsRatio = rmsList.back() == 0.f ? 1.f : rms / rmsList.back();
        rmsList.push_back(rms);

        if (rmsRatio > factor && rms > threshold) {
            transientFound = true;
            rmsWindow = tempRms; // Save the window that triggered
            tInitSample = samp;
            break;
        }
    }

    if (!transientFound) {
        // No transient, just find the first peak in a small window
        int peakSample = Sitrano::findPeakSample(aud, 0, std::min((int)aud.size(), 4096), true);
        range.endSample = Sitrano::findPreviousZero(aud, peakSample);
        return range;
    }

    int peakIndexInWindow = Sitrano::findAbsPeakIndex(rmsWindow);
    int peakSample = tInitSample + peakIndexInWindow;
    float peakValue = std::abs(aud[peakSample]);

    // Find the *actual* number of samples in that first window
    int foundWindowNumSamples = 0;
    for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + tInitSample) < aud.size(); rmsSamp++) {
        foundWindowNumSamples++;
    }

    // Check if the peak is the last sample of that window
    if (foundWindowNumSamples > 0 && peakIndexInWindow == (foundWindowNumSamples - 1))
    {
        // The peak is at the end. We must search subsequent windows.
        for (int samp = tInitSample + rmsHopLength; samp < aud.size(); samp += rmsHopLength)
        {
            std::vector<float> nextRmsWindow(rmsSize, 0.f);
            int currentWindowNumSamples = 0;
            for (int rmsSamp = 0; rmsSamp < rmsSize && (rmsSamp + samp) < aud.size(); rmsSamp++)
            {
                nextRmsWindow[rmsSamp] = aud[samp + rmsSamp];
                currentWindowNumSamples++;
            }

            if (currentWindowNumSamples == 0) {
                break; // Reached the end of the audio. The last peak we found stands.
            }

            int nextPeakIndexInWindow = Sitrano::findAbsPeakIndex(nextRmsWindow);

            if (nextPeakIndexInWindow >= currentWindowNumSamples) {
                continue;
            }

            float nextPeakValue = std::abs(nextRmsWindow[nextPeakIndexInWindow]);

            if (nextPeakValue > peakValue)
            {
                peakValue = nextPeakValue;
                peakSample = samp + nextPeakIndexInWindow;
                if (nextPeakIndexInWindow != (currentWindowNumSamples - 1)) break;
            }
            else break;
        }
    }

    range.initSample = Sitrano::findPreviousZero(aud, peakSample);
    range.endSample = Sitrano::findNextZero(aud, peakSample);
    return range;
}