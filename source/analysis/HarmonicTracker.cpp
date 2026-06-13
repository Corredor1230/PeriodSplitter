#include"HarmonicTracker.h"
#include"support/AudioChunkStrategy.h"
#include"support/PeriodLoopStrategy.h"
#include"support/SinglePeriodStrategy.h"

HarmonicTracker::HarmonicTracker(
    const Sihat::AnalysisUnit& a,
    const Sihat::AnalysisConfig& conf,
    const std::vector<Sihat::Peak>& top,
    const std::vector<uint32_t>& sampleList,
    const int startSample, const float pitch) :
    unit(a),
    start(startSample),
    sr(unit.sampleRate),
    config(conf),
    sList(sampleList),
    nfft(conf.nfft),
    settings(conf.hSettings),
    tFreqs(top),
    f0(pitch)
{
    window.resize(nfft);
    checker.resize(nfft * 2);
    for (int i = 0; i < nfft; i++)
    {
        window[i] = 0.0;
    }
    initFFTW();

    switch (settings.style)
    {
    case Sihat::WindowStyle::periodLoop:
        mWindowStrategy = std::make_unique<PeriodLoopStrategy>();
        break;
    case Sihat::WindowStyle::singlePeriod:
        mWindowStrategy = std::make_unique<SinglePeriodStrategy>();
        break;
    case Sihat::WindowStyle::audioChunk:
    default:
        mWindowStrategy = std::make_unique<AudioChunkStrategy>();
        break;
    }
}

void HarmonicTracker::initFFTW() {
    input = (float*)fftwf_malloc(sizeof(float) * nfft);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(nfft, input, output, FFTW_MEASURE);
}

void HarmonicTracker::applyHann(float* data, int size) {
    for (int i = 0; i < size; ++i) {
        float w = 0.5f * (1.0f - std::cos(Sihat::TWO_PI * i / (size - 1)));
        data[i] *= w;
    }
}

double HarmonicTracker::interpolatePeak(int k, const std::vector<double>& mags)
{
    if (k <= 0 || k >= (int)mags.size() - 1) return 0.0;
    double m1 = mags[k - 1], m0 = mags[k], p1 = mags[k + 1];
    double denom = (m1 - 2 * m0 + p1);
    if (fabs(denom) < 1e-12) return 0.0;
    return 0.5 * (m1 - p1) / denom;
}

Sihat::HarmonicResults HarmonicTracker::analyze()
{
    LoopParameters params = mWindowStrategy->getLoopParameters(sList, unit, config, start);

    hResults.amps.resize(config.numHarmonics);
    hResults.freqs.resize(config.numHarmonics);
    hResults.phases.resize(config.numHarmonics);

    size_t numFrames = params.numFrames;
    int frameStep = params.frameStep;

    float ampThresh = Sihat::dbToAmp(-50.0);

    int peakSample = findPeakSample();

    for (size_t i = 0; i < numFrames - (frameStep + 1); i += frameStep) {

        int periodLength{ 0 };

        bool success = mWindowStrategy->processFrame(
            i, frameStep, config, unit, sList, hResults, input, periodLength, start
        );

        if (!success) continue;
        if (settings.applyHanning) applyHann(input, nfft);

        fftwf_execute(plan);

        for (int h = 1; h <= tFreqs.size(); ++h) {
            if (tFreqs[h - 1].freq <= 20.0) continue;
            float targetFreq = tFreqs[h - 1].freq;
            Sihat::BinFreq binFreq = settings.tolInHz ? Sihat::findPeakWithinTolerance(targetFreq, settings.tolInHz, nfft,
                sr, output, true) : Sihat::findPeakWithinTolerance(targetFreq, settings.toleranceValue, nfft, sr, output, false);

            if (binFreq.bin < nfft / 2 + 1) {
                Sihat::FreqUnit fUnit{ binFreq, 0.0, 0.0, 0.0 };
                //int binRange = 5; //Gotta improve this eventually, make algorithm

                fUnit = Sihat::findPeak(fUnit.bin, output, nfft, unit.sampleRate);

                if (std::isnan(fUnit.amp)) 
                break;

                for (int step = 0; step < frameStep; step++)
                {
                    if (fUnit.amp > 1.0)
                        std::cout << "i: " << i << " , h: " << h << " , val: " << fUnit.amp << '\n';

                    float duringEnv = 1.0;
                    if (i < hResults.finalSamples.size()){
                        if (hResults.finalSamples[i] - start < peakSample && peakSample != 0)
                        {
                            float rate = static_cast<float>(hResults.finalSamples[i] - start) / static_cast<float>(peakSample); 
                            duringEnv = Sihat::normLogistic(rate, 2.0);
                        }
                    }

                    hResults.amps[h - 1].push_back(Sihat::ampToDb(fUnit.amp * duringEnv));
                    hResults.phases[h - 1].push_back(fUnit.pha);
                    hResults.freqs[h - 1].push_back(fUnit.bin.freq / f0);
                }
            }
        }
    }
    if (hResults.finalSamples.size() > hResults.amps[0].size())
    {
        hResults.finalSamples.resize(hResults.amps[0].size());
    }
    else if (hResults.finalSamples.size() < hResults.amps[0].size())
    {
        for (int i = hResults.finalSamples.size(); i < hResults.amps[0].size(); i++)
        {
            std::cerr << "This is wrong\n";
        }
    }

    //Extract RMS of first second of audio
    float sum = 0.0;
    for (int i = start; i < start + sr; ++i)
    {
        sum += unit.soundFile[i] * unit.soundFile[i];
    }
    hResults.rms = std::sqrt(sum / sr);
    return hResults;
}

int HarmonicTracker::findPeakSample()
{
    int rmsSize = static_cast<int>(sr / 200.0); // 10ms window
    int hopSize = rmsSize;
    int maxScanDuration = static_cast<int>(sr); // 1 second of audio
    
    // Determine the absolute hard limit of the audio file to prevent crashes
    int totalSamples = static_cast<int>(unit.soundFile.size());
    int safeEnd = std::min(start + maxScanDuration, totalSamples - rmsSize);

    float maxRMS = -1.0f;
    int bestSampleIndex = start;
    float rmsThresh = 0.75;

    std::vector<float> rmsList;

    // Scan through the file in safe hops
    for (int frameStart = start; frameStart <= safeEnd; frameStart += hopSize)
    {
        float sum = 0.0f;
        for (int i = 0; i < rmsSize; i++)
        {
            float sample = unit.soundFile[frameStart + i];
            sum += sample * sample; 
        }
        
        float currentRMS = std::sqrt(sum / static_cast<float>(rmsSize));
        rmsList.push_back(currentRMS);
        
        if (currentRMS > maxRMS)
        {
            maxRMS = currentRMS;
            //bestSampleIndex = frameStart;
        }
    }

    for (int i = 0; i < rmsList.size(); i++){
        if (rmsList[i] >= rmsThresh* maxRMS){
            bestSampleIndex = i * hopSize + start; 
            break;
        }
    }

    return bestSampleIndex;
}