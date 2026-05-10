#include"TransientRecognition.h"
#include"dsp/ZPFilter.h"
#include<deque>
#include<numeric>
#include"dsp/MorletWavelet.h"

Transient::Transient(
    const Sihat::AnalysisUnit& unit,
    const Sihat::TransientSettings& conf,
    const Sihat::TransientFFTSettings& ffts,
    const float pitch
) :
    tSettings(conf),
    pitch(pitch),
    sampleRate(unit.sampleRate),
    aud(unit.soundFile),
    factor(conf.transientFactor),
    threshold(conf.transientThreshold),
    corrThreshold(conf.tCorrelationThreshold),
    preAttack(conf.preAttack),
    sr(unit.sampleRate),
    nfft(ffts.nfft),
    hop(ffts.hopSize),
    useFFT(ffts.useFFT),
    flatnessThresh(ffts.flatnessThreshold),
    correlationOffset(2500)
{

    plan, input, output = nullptr;

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

    initFFTW();
}

Transient::~Transient()
{
    if (plan != nullptr) {
        fftwf_destroy_plan(plan);
    }

    if (input != nullptr) {
        fftwf_free(input); 
    }

    if (output != nullptr) {
        fftwf_free(output);
    }
}

void Transient::initFFTW()
{
    input = (float*)fftwf_malloc(sizeof(float) * nfft);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(nfft, input, output, FFTW_MEASURE);
}

Sihat::TransientResults Transient::findStartTransient()
{
    Sihat::TransientResults results;
    Sihat::SampleRange range{0, 0};
    if (!useFFT)
    {
        range = findFromRMS();
        bool transientFailed = false;
        if (range.initSample == 0 || range.endSample == 0)
        { 
            transientFailed = true;
            range = findFromFFT();
        }
    }
    else
    {
        range = findFromFFT();
    }

    range = findWithCrossCorrelation(correlationOffset, range.initSample);

    std::vector<float> transient(range.endSample - range.initSample);

    for (int i = 0; i < transient.size(); ++i)
    {
        transient[i] = aud[range.initSample + i];
    }

    MorletWavelet wavelet(transient, sr);

    std::vector<float> wTransform(transient.size());

    for (int i = 0; i < wTransform.size(); ++i)
    {
        float midiLow = Sihat::freqToMidi(pitch / 2.f);
        float midiHigh = Sihat::freqToMidi(20000.f);

        float midiRange = midiHigh - midiLow;
        float midiStep = midiRange / static_cast<float>(transient.size());

        wTransform[i] = Sihat::midiToFreq(midiLow + midiStep * i);
    }

    std::vector<Sihat::VariableRatePartial> scalogram = wavelet.analyze(64, 2048, pitch * 1.0, 20000.f);

    double sourceSumSq = 0.0;
    int rangeLen = range.endSample - range.initSample;
    if (rangeLen > 0) {
        for (int i = 0; i < rangeLen; ++i) {
            float sample = aud[range.initSample + i];
            sourceSumSq += sample * sample;
        }
    }
    double sourceRMS = (rangeLen > 0) ? std::sqrt(sourceSumSq / rangeLen) : 0.0;

    // 2. Calculate the estimated RMS of your SCALOGRAM (The Model)
    // We assume the partials are mostly incoherent (energies sum up).
    // Power of a sine wave = (Amplitude^2) / 2.
    
    double modelTotalPower = 0.0;
    
    for (const auto& partial : scalogram) {
        if (partial.data.empty()) continue;

        double partialSumSq = 0.0;
        for (float val : partial.data) {
            partialSumSq += val * val;
        }
        
        // Average squared amplitude for this partial
        double partialMeanSq = partialSumSq / partial.data.size();
        
        // Convert Peak Amplitude to Power (RMS^2) -> A^2 / 2
        double partialPower = partialMeanSq / 2.0; 
        
        modelTotalPower += partialPower;
    }
    
    double modelRMS = std::sqrt(modelTotalPower);

    // 3. Determine the Gain Factor
    // Prevent division by zero
    float gain = 1.0f;
    if (modelRMS > 0.00001) {
        gain = static_cast<float>(sourceRMS / modelRMS);
    }
    
    // 4. Apply Gain to the Scalogram
    // This forces the "Synthesized Transient" to have the exact same energy as the "Recorded Transient"
    for (auto& partial : scalogram) {
        for (float& val : partial.data) {
            val *= gain;
        }
    }

    std::cout << "Transient Calibration | Source RMS: " << sourceRMS 
            << " | Model RMS: " << modelRMS 
            << " | Applied Gain: " << gain << std::endl;

    results.range = range;
    results.scalogram = scalogram;

    //Extract original audio transient RMS for resynthesis
    float sum = 0.0;
    for (int i = results.range.initSample; i < results.range.endSample; ++i)
    {
        sum += aud[i] * aud[i];
    }
    results.rms = std::sqrt(sum / (results.range.endSample - results.range.initSample));

    return results;
}

Sihat::SampleRange Transient::findFromRMS()
{
    Sihat::SampleRange range;
    range.initSample = tSettings.tStartSample;

    std::vector<float> rmsList;
    bool transientFound = false;
    int tInitSample = 0;

    std::vector<float> rmsWindow(rmsSize); // This will hold the transient

    for (int samp = tSettings.tStartSample; samp < (aud.size() / 4); samp += rmsHopLength)
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
            rmsWindow = tempRms;
            tInitSample = samp;
            break;
        }
    }

    if (!transientFound) {
        // No transient, just find the first peak in a small window
        int peakSample = Sihat::findPeakSample(aud, 0, std::min((int)aud.size(), 8192), true);
        range.endSample = Sihat::findPreviousZero(aud, peakSample);
        return range;
    }

    range.initSample = Sihat::findPreviousZero(aud, tInitSample);
    range.endSample = Sihat::findNextZero(aud, tInitSample + rmsSize);
    return range;
}

// Rewritten transient detector using spectral flatness
Sihat::SampleRange Transient::findFromFFT()
{
    Sihat::SampleRange range{0, 0};

    // 1. Basic Safety Checks
    if (nfft <= 2 || hop <= 0 || aud.empty()) return range;

    // Define search window
    // If RMS failed, start from 0. Otherwise, look a bit before the RMS detection.
    int startSamp = 0;
    //int startSamp = transientFailed ? 0 : std::max(0, rmsRange.initSample - (4 * hop));
    const int lastStart = static_cast<int>(aud.size()) - nfft;

    if (startSamp > lastStart) return range;

    // 2. Prepare Windowing (Hanning)
    std::vector<float> window(nfft);
    for (int i = 0; i < nfft; ++i) {
        window[i] = 0.5f * (1.0f - std::cos(2.0f * static_cast<float>(Sihat::PI) * i / (nfft - 1)));
    }

    // 3. Setup Frequency Limits
    // Ignore DC/Rumble (< 40Hz) and Super-high noise (> 12kHz)
    const float minFreq = 40.0f;
    const float maxFreq = 12000.0f; 
    const int firstBin = std::max(1, static_cast<int>(minFreq * nfft / sr));
    const int lastBin  = std::min(nfft / 2 - 1, static_cast<int>(maxFreq * nfft / sr));

    // 4. Flux State Variables
    // Stores magnitude of previous frame for subtraction
    std::vector<float> prevMag(nfft / 2 + 1, 0.0f); 
    
    // Adaptive History (Rolling Average)
    std::deque<float> fluxHistory;
    const int historySize = 8; // Look at last ~8 frames for average
    
    // Parameters
    const float sensitivity = 2.0f; // Multiplier: Current must be 2x the average
    const float noiseFloor = 0.05f; // Absolute minimum flux to trigger (normalized)

    bool foundTransient = false;
    int transientFrameStart = -1;

    // 5. Analysis Loop
    for (int i = startSamp; i <= lastStart; i += hop)
    {
        // A. Load and Window Audio
        for (int j = 0; j < nfft; ++j) {
            input[j] = aud[i + j] * window[j];
        }

        fftwf_execute(plan);

        // B. Calculate Spectral Flux
        float currentFlux = 0.0f;

        for (int b = firstBin; b <= lastBin; ++b)
        {
            float real = output[b][0];
            float imag = output[b][1];
            
            // Calculate Magnitude
            float mag = std::sqrt(real * real + imag * imag);

            // Rectified Difference: We only care if energy INCREASED
            float diff = mag - prevMag[b];
            if (diff > 0) {
                currentFlux += diff;
            }

            // Update history for next frame
            prevMag[b] = mag;
        }

        // C. Normalize Flux (Crucial for consistent thresholds)
        currentFlux /= static_cast<float>(nfft);

        // D. Adaptive Threshold Logic
        float localAverage = 0.0f;
        if (!fluxHistory.empty()) {
            localAverage = std::accumulate(fluxHistory.begin(), fluxHistory.end(), 0.0f) 
                           / static_cast<float>(fluxHistory.size());
        }

        // The dynamic threshold: Relative increase + Absolute floor
        float threshold = (localAverage * sensitivity) + noiseFloor;

        // Check Trigger
        // (fluxHistory.size() check ensures we don't trigger on the very first frame)
        if (currentFlux > threshold && fluxHistory.size() >= 3)
        {
            float windowsum = 0.0;
            for (int r = 0; r < nfft; r++)
            {
                windowsum += std::pow(aud[i + r], 2);
            }
            float windowrms = std::sqrt(windowsum / (float)nfft);

            if (windowrms > 0.2)
            {
                foundTransient = true;
                transientFrameStart = findFirstAboveThreshold(i, 0.3);
                break; // Stop at first strong onset
            }
        }

        // E. Update History
        fluxHistory.push_back(currentFlux);
        if (fluxHistory.size() > historySize) {
            fluxHistory.pop_front();
        }
    }

    // 6. Return Result
    if (foundTransient)
    {
        // Pad the result slightly to ensure we capture the very start of the rise
        range.initSample = Sihat::findNearestZero(aud, std::max(0, transientFrameStart - hop));
        range.endSample = Sihat::findNearestZero(aud, std::min(static_cast<int>(aud.size()), range.initSample + nfft));
    }
    else
    {
        // Fallback to RMS if FFT analysis was inconclusive
        range = {0, 0};
    }

    return range;
}

int Transient::findFirstAboveThreshold(int startSample, float thresh)
{
    int consecutiveSamples = 0;
    int firstAbove = 0;
    for (int i = startSample; i < aud.size(); ++i)
    {
        if (aud[i] > thresh) ++consecutiveSamples;
        else consecutiveSamples = 0;

        if (consecutiveSamples > 3)
        {
            firstAbove = i - 2;
            break;
        }
    }    
    
    if (firstAbove > aud.size() / 2) return 0;
    else return firstAbove;
}

Sihat::SampleRange Transient::findWithCrossCorrelation(int offset, int firstSample)
{
    int withOffset = firstSample + offset;
    int expectedPeriod = static_cast<int>(sr / pitch);
    int slidingSize = expectedPeriod - static_cast<int>(expectedPeriod * 0.1);
    float maxPitchRatio = 8.0;

    ZPFilter filter;

    int corrStart = Sihat::findPreviousZero(aud, Sihat::findPeakSample(aud, withOffset, withOffset + expectedPeriod * 2, true));
    int firstWindowEnd = Sihat::findNearestZero(aud, corrStart + expectedPeriod);

    std::vector<float> window(firstWindowEnd - corrStart);

    //We copy data from audio to the window
    for (int i = 0; i < window.size(); ++i)
    {
        window[i] = aud[corrStart + i];
    }

    float squareA = 0.f;

    window = filter.filtfilt(window, sr, pitch * maxPitchRatio);

    for (float w : window) squareA += w * w;


    //Yes, we are going to have two different vectors that copy data from audio. Yes, this is costly, but it's best if we want to filter out higher frequency content.
    std::vector<float> comparative(window.size());

    bool correlationFound = false;
    std::deque<float> corrList;
    int earliestCorrelation = 0;
    for (int w = corrStart; w > firstSample; --w)
    {
        int compStart = w - slidingSize;
        for(int i = 0; i < comparative.size(); ++i)
        {
            comparative[i] = aud[compStart + i];
        }
        comparative = filter.filtfilt(comparative, sr, pitch * maxPitchRatio);

        float correlationValue = 0.f;

        float numerator = 0.f;
        float squareB = 0.f;

        for (int i = 0; i < window.size(); ++i)
        {
            float a = window[i];
            float b = comparative[i];
            numerator += a * b;
            squareB += b * b;
        }

        float denominator = std::sqrt(squareA * squareB);

        correlationValue = (denominator == 0.f) ? 0.f : numerator / denominator;
        bool aboveThreshold = false;

        corrList.emplace_front(correlationValue);

        if (corrList.size() < 3)
        {
            continue;
        }

        int lc = corrList.size() - 1;

        aboveThreshold = corrList[1] > corrThreshold;
        bool isPeak = (corrList[1] > corrList[0]) && (corrList[1] > corrList[2]);

        if (isPeak && aboveThreshold)
        {
            earliestCorrelation = Sihat::findNearestZero(aud, compStart);
            for (int i = 0; i < window.size(); ++i)
            {
                window[i] = aud[earliestCorrelation + i];
            }
            filter.filtfilt(window, sr, pitch * maxPitchRatio);
            squareA = 0.f;
            for (float k : window) squareA += k * k;

            w -= slidingSize;
        }

    }

    Sihat::SampleRange out{firstSample, earliestCorrelation + 1};

    //Maybe we can consider a Haar Wavelet transform for this.
    return out;
}