#include"HarmonicTracker.h"
#include"support/AudioChunkStrategy.h"
#include"support/PeriodLoopStrategy.h"
#include"support/SinglePeriodStrategy.h"

HarmonicTracker::HarmonicTracker(
    const Sitrano::AnalysisUnit& a,
    const Sitrano::AnalysisConfig& conf,
    const std::vector<Sitrano::Peak>& top,
    const std::vector<uint32_t>& sampleList) :
    unit(a),
    config(conf),
    sList(sampleList),
    settings(conf.hConfig),
    tFreqs(top)
{
    nfft = config.nfft;
    window.resize(nfft);
    checker.resize(nfft * 2);
    for (int i = 0; i < nfft; i++)
    {
        window[i] = 0.0;
    }
    initFFTW();

    switch (settings.style)
    {
    case Sitrano::WindowStyle::periodLoop:
        mWindowStrategy = std::make_unique<PeriodLoopStrategy>();
        break;
    case Sitrano::WindowStyle::singlePeriod:
        mWindowStrategy = std::make_unique<SinglePeriodStrategy>();
        break;
    case Sitrano::WindowStyle::audioChunk:
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
        float w = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (size - 1)));
        data[i] *= w;
    }
}

//void HarmonicTracker::analyze() {
//    LoopParameters params = mWindowStrategy->getLoopParameters(results, unit);
//    size_t numFrames = params.numFrames;
//    int frameStep = params.frameStep;
//
//    float ampThresh = Sitrano::dbToAmp(-50.0);
//
//    for (size_t i = 0; i < numFrames - (frameStep + 1); i+=frameStep) {
//
//        int periodLength{ 0 };
//
//        bool success = mWindowStrategy->processFrame(
//            i, frameStep, unit, results, input, periodLength
//        );
//
//        if (!success) continue;
//        if (settings.applyHanning) applyHann(input, unit.nfft);
//
//        fftwf_execute(plan);
//        
//        for (int h = 1; h <= tFreqs.size(); ++h) {
//            if (tFreqs[h - 1].freq <= 20.0) continue;
//            float targetFreq = tFreqs[h - 1].freq;
//            float binFreq = 0.f;
//            binFreq = Sitrano::findPeakWithinTolerance(targetFreq, settings.toleranceValue, unit.nfft, unit.sampleRate, output);
//            
//            int bin = static_cast<int>(std::round(targetFreq * unit.nfft / float(unit.sampleRate)));
//
//            if (bin < unit.nfft / 2 + 1) {
//                float real = output[bin][0];
//                float imag = output[bin][1];
//
//                float mag = std::sqrt(real * real + imag * imag);
//                float amp = Sitrano::mag_to_amp(mag, unit.nfft);
//                float phase = std::atan2(imag, real);  // radians, range [-?, ?]
//                
//                for (int step = 0; step < frameStep; step++)
//                {
//                    results.amps[h - 1].push_back(amp);
//                    results.phases[h - 1].push_back(phase);
//                    results.freqs[h - 1].push_back(binFreq);
//                    
//                }
//            }
//        }
//    }
//    if (results.finalSamples.size() > results.amps[0].size())
//    {
//        results.finalSamples.resize(results.amps[0].size());
//    }
//    else if (results.finalSamples.size() < results.amps[0].size())
//    {
//        for (int i = results.finalSamples.size(); i < results.amps[0].size(); i++)
//        {
//            std::cerr << "This is wrong\n";
//        }
//    }
//}

double HarmonicTracker::interpolatePeak(int k, const std::vector<double>& mags)
{
    if (k <= 0 || k >= (int)mags.size() - 1) return 0.0;
    double m1 = mags[k - 1], m0 = mags[k], p1 = mags[k + 1];
    double denom = (m1 - 2 * m0 + p1);
    if (fabs(denom) < 1e-12) return 0.0;
    return 0.5 * (m1 - p1) / denom;
}

Sitrano::HarmonicResults HarmonicTracker::getEnvelopes()
{
    LoopParameters params = mWindowStrategy->getLoopParameters(sList, unit, config);
    //Sitrano::HarmonicResults hResults;

    hResults.amps.resize(config.numHarmonics);
    hResults.freqs.resize(config.numHarmonics);
    hResults.phases.resize(config.numHarmonics);

    size_t numFrames = params.numFrames;
    int frameStep = params.frameStep;

    float ampThresh = Sitrano::dbToAmp(-50.0);

    for (size_t i = 0; i < numFrames - (frameStep + 1); i += frameStep) {

        int periodLength{ 0 };

        bool success = mWindowStrategy->processFrame(
            i, frameStep, config, unit, sList, hResults, input, periodLength
        );

        if (!success) continue;
        if (settings.applyHanning) applyHann(input, nfft);

        fftwf_execute(plan);

        for (int h = 1; h <= tFreqs.size(); ++h) {
            if (tFreqs[h - 1].freq <= 20.0) continue;
            float targetFreq = tFreqs[h - 1].freq;
            float binFreq = 0.f;
            binFreq = Sitrano::findPeakWithinTolerance(targetFreq, settings.toleranceValue, nfft, unit.sampleRate, output);

            int bin = static_cast<int>(std::round(targetFreq * nfft / float(unit.sampleRate)));

            if (bin < nfft / 2 + 1) {
                float real = output[bin][0];
                float imag = output[bin][1];

                float mag = std::sqrt(real * real + imag * imag);
                float amp = Sitrano::mag_to_amp(mag, nfft);
                float phase = std::atan2(imag, real);  // radians, range [-?, ?]

                for (int step = 0; step < frameStep; step++)
                {
                    hResults.amps[h - 1].push_back(amp);
                    hResults.phases[h - 1].push_back(phase);
                    hResults.freqs[h - 1].push_back(binFreq);

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

    return hResults;
}