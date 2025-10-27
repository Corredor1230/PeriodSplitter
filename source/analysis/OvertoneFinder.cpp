#include "OvertoneFinder.h"

OvertoneFinder::OvertoneFinder(const Sitrano::AnalysisUnit& unit,
    const Sitrano::AnalysisConfig& conf) : 
    unit(unit),
    config(conf),
    settings(conf.oConfig),
    N(config.nfft * 4),
    Nout(N / 2 + 1)
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

    // --- 1. PREPARE THE INPUT BUFFER ---
    // We use the pre-allocated 'm_fft_in' member buffer

    // Find how much data to copy: either the whole signal or just enough to fill nfft
    size_t numToCopy = std::min(checkSignal.size(), (size_t)N);

    // Clear the buffer (in case previous run left data)
    std::fill_n(input, N, 0.0f);

    // Copy the signal. If checkSignal.size() < m_nfft, this is zero-padded.
    std::copy(checkSignal.begin(), checkSignal.begin() + numToCopy, input);

    // Apply the Hanning window to the entire nfft buffer
    // (Using the safer vector-based function from our last conversation)
    Sitrano::applyHann<float>(input, N);

    // --- 2. EXECUTE FFT ---
    // We re-use the *exact same plan* every time. This is much faster.
    fftwf_execute(plan);

    // --- 3. CALCULATE MAGNITUDES ---
    std::vector<double> mags(Nout);
    for (int k = 0; k < Nout; ++k) {
        double re = output[k][0];
        double im = output[k][1];
        mags[k] = std::sqrt(re * re + im * im);
    }

    // --- 4. FIND & INTERPOLATE PEAKS ---
    std::vector<Sitrano::Peak> peaks;
    peaks.reserve(Nout); // Good practice
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

    // --- 5. CLUSTER PEAKS TO FIND HARMONICS ---

    // Implement the tolerance logic correctly
    auto withinTolerance = [&](double baseFreq, double compFreq) {
        if (settings.useTolerance) {
            double tolHz = Sitrano::cents_to_hz(baseFreq, settings.toleranceValue);
            return std::abs(compFreq - baseFreq) <= tolHz;
        }
        else {
            // Assume 'tolerance' is in Hz if 'useCentTolerance' is false
            return std::abs(compFreq - baseFreq) <= settings.toleranceValue;
        }
        };

    std::vector<Sitrano::Peak> merged;
    merged.reserve(outHarm);

    // This is a much safer way to track processed peaks
    std::vector<bool> peak_processed(peaks.size(), false);

    for (int i = 0; i < peaks.size() && merged.size() < outHarm; ++i) {
        if (peak_processed[i]) {
            continue; // This peak was already added to a cluster
        }

        // This is our new "seed" peak
        const auto& seedPeak = peaks[i];
        peak_processed[i] = true;

        // Apply filters to the seed peak
        if (seedPeak.freq < pitch * 0.9 ||
            seedPeak.freq > 20000.f ||
            seedPeak.amp < 1.0e-6) {
            continue; // This peak is invalid, skip it
        }

        Sitrano::Peak cluster = seedPeak;

        // Implement 'sumAmplitudesInCluster' logic
        if (settings.sumAmplitudes) {
            // Iterate over *remaining* peaks to find matches
            for (int j = i + 1; j < peaks.size(); ++j) {
                if (peak_processed[j]) continue;

                if (withinTolerance(cluster.freq, peaks[j].freq)) {
                    // It's in the cluster. Sum its amplitude.
                    cluster.amp += peaks[j].amp;
                    cluster.mag += peaks[j].mag; // Also sum mag
                    peak_processed[j] = true; // Mark as processed
                }
            }
        }

        merged.push_back(cluster);
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