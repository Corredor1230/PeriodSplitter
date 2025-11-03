#include "OvertoneFinder.h"
#include<algorithm>

OvertoneFinder::OvertoneFinder(const Sitrano::AnalysisUnit& unit,
    const Sitrano::AnalysisConfig& conf) : 
    unit(unit),
    config(conf),
    settings(conf.oConfig),
    N(config.oConfig.fftSize),
    Nout(N / 2 + 1),
    absThreshold(conf.oConfig.setAbsoluteThreshold),
    threshold(conf.oConfig.overtoneThreshold)
{
    initFFTW();
}

void OvertoneFinder::initFFTW() {
    input = (float*)fftwf_malloc(sizeof(float) * config.oConfig.fftSize);
    output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (config.oConfig.fftSize / 2 + 1));
    plan = fftwf_plan_dft_r2c_1d(config.oConfig.fftSize, input, output, FFTW_MEASURE);
}

// This is now a member function of OvertoneFinder
std::vector<Sitrano::Peak> OvertoneFinder::getRelevantOvertones(
    const std::vector<double>& checkSignal,
    float pitch)
{
    const int outHarm = config.numHarmonics;
    N = config.oConfig.fftSize;
    const int Nout = N / 2 + 1;
    size_t numToCopy = std::min(checkSignal.size(), (size_t)N);

    std::fill_n(input, N, 0.0f);
    std::copy(checkSignal.begin(), checkSignal.begin() + numToCopy, input);

    Sitrano::applyHann<float>(input, N);

    fftwf_execute(plan);

    // CALCULATE MAGNITUDES
    std::vector<double> mags(Nout);
    for (int k = 0; k < Nout; ++k) {
        double re = output[k][0];
        double im = output[k][1];
        mags[k] = std::sqrt(re * re + im * im);
    }

    // FIND & INTERPOLATE PEAKS 
    std::vector<Sitrano::Peak> peaks;
    peaks.reserve(Nout); 
    for (int k = 1; k < Nout - 1; ++k) {
        // Find local maxima (basic peak-picking)
        if (mags[k] > mags[k - 1] && mags[k] > mags[k + 1])
        {
            // Note: Your Sitrano::interp_delta probably does quadratic interpolation,
            // which is a great approach.
            double delta = Sitrano::interp_delta(k, mags);
            double f = (k + delta) * unit.sampleRate / double(N);
            double amp = Sitrano::mag_to_amp(mags[k], N);
            peaks.push_back({ f, mags[k], amp });
        }
    }

    // Sort all found peaks by amplitude (loudest first)
    std::sort(peaks.begin(), peaks.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
        return a.amp > b.amp;
        });

    std::vector<Sitrano::Peak> merged;
    merged.reserve(outHarm);

    // This is a much safer way to track processed peaks
    std::vector<bool> peak_processed(peaks.size(), false);

    float ampThresh = 0.5;
    if (absThreshold)
    {
        ampThresh = Sitrano::dbToAmp(threshold);
    }
    else
    {
        float peakDb = Sitrano::ampToDb(peaks[0].amp);
        ampThresh = Sitrano::dbToAmp(peakDb - std::abs(threshold));
    }

    for (int i = 0; i < peaks.size() && merged.size() < outHarm; ++i) {
        if (peak_processed[i]) {
            continue; // This peak was already added to a cluster
        }

        // This is our new "seed" peak
        const auto& seedPeak = peaks[i];
        float tolInHz = Sitrano::cents_to_hz(seedPeak.freq, settings.toleranceValue);
        peak_processed[i] = true;
        bool withinTolerance = false;

        // Apply filters to the seed peak
        if (seedPeak.freq < pitch * 0.9 ||
            seedPeak.freq > 20000.f ||
            seedPeak.amp < ampThresh) {
            continue; // This peak is invalid, skip it
        }

        Sitrano::Peak cluster = seedPeak;

        for (int i = 0; i < merged.size(); i++)
        {
            if (Sitrano::withinTolerance(merged[i].freq, seedPeak.freq, tolInHz))
            {
                withinTolerance = true;
            }
        }

        if (!withinTolerance) merged.push_back(cluster);
    }

    // --- 6. FINAL SORTING & CLEANUP ---

    // If we summed amplitudes, the list is no longer sorted by amp.
    // We must re-sort to get the Top N.
    if (settings.sumAmplitudes) {
        std::sort(merged.begin(), merged.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
            return a.amp > b.amp;
            });
    }

    // Now, resize to the desired number of harmonics
    if (merged.size() > outHarm) {
        merged.resize(outHarm);
    }

    // Finally, sort by frequency for a clean output
    std::sort(merged.begin(), merged.end(), [](const Sitrano::Peak& a, const Sitrano::Peak& b) {
        return a.freq < b.freq;
        });

    // No need to destroy the plan, it lives in the class instance
    return merged;
}