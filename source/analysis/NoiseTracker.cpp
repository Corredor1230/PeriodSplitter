#include "NoiseTracker.h"

NoiseTracker::NoiseTracker(const Sitrano::AnalysisUnit& ana, 
	Sitrano::Results& r,
	float startFreq):
	unit(ana),
	r(r),
	startFreq(startFreq)
{
	float f0 = startFreq;
	while (f0 < ana.sampleRate / 2.0 || f0 < 20000.0) {
		float f_low = f0 / pow(2.0, 1.0 / 6.0);
		float f_high = f0 * pow(2.0, 1.0 / 6.0);
		bands.push_back({ f_low, f_high, f0 });
		f0 *= pow(2.0, 1.0 / 3.0);
	}

	// FFTW Initialization
	N = unit.nfft / 4;
	fft_in = fftwf_alloc_real(N);
	fft_out = fftwf_alloc_complex(N / 2 + 1);

	plan = fftwf_plan_dft_r2c_1d(N, fft_in, fft_out, FFTW_MEASURE);
}

NoiseTracker::~NoiseTracker()
{
	fftwf_destroy_plan(plan);
	fftwf_free(fft_in);
	fftwf_free(fft_out);
}

void NoiseTracker::applyFrameTable(std::vector<int> table)
{
	frameTable = table;
	useFrameTable = true;
}

void NoiseTracker::analyze() {
    // --- Initialization (mostly unchanged) ---
    const int N = unit.nfft;
    const int hop = unit.hopSize;
    const std::vector<float>& signal = unit.soundFile;
    const float sampleRate = unit.sampleRate;

    int numFrames = useFrameTable ? frameTable.size() : (signal.size() - N) / hop;
    if (numFrames <= 0) return;
    const float df = sampleRate / N;

    // --- 1. Restructure Data Format ---
    // The outer vector now represents BANDS, the inner vector represents FRAMES (time).
    r.noise.assign(bands.size(), std::vector<float>(numFrames, 0.0f));

    // --- Main Analysis Loop (calculates power) ---
    for (int i = 0; i < numFrames; ++i) {
        int startSample = useFrameTable ? frameTable[i] : i * hop;
        if (startSample + N > signal.size()) {
            // Trim all band vectors if the signal ends unexpectedly
            for (auto& band_vec : r.noise) {
                band_vec.resize(i);
            }
            break;
        }

        // Get frame, apply window, and run FFT (unchanged)
        std::vector<float> frame(signal.begin() + startSample, signal.begin() + startSample + N);
        Sitrano::applyWindow(frame);
        std::copy(frame.begin(), frame.end(), fft_in);
        fftwf_execute(plan);

        // Sum power into the correct bands
        for (int k = 0; k < N / 2 + 1; ++k) {
            float freq = k * df;
            float real = fft_out[k][0];
            float imag = fft_out[k][1];
            float power = real * real + imag * imag;

            for (int b = 0; b < bands.size(); ++b) {
                if (freq >= bands[b].f_low && freq < bands[b].f_high) {
                    // Assign to transposed matrix: results[band][frame]
                    r.noise[b][i] += power;
                    break;
                }
            }
        }
    }

    // --- 2. Normalization Step ---
    // This runs after all power values have been calculated for all frames.
    for (auto& band_envelope : r.noise) {
        // Find the peak energy in this band's envelope using std::max_element
        auto max_iterator = std::max_element(band_envelope.begin(), band_envelope.end());

        // Proceed only if the band is not empty
        if (max_iterator != band_envelope.end()) {
            float maxEnergy = *max_iterator;

            // Normalize the envelope if the peak is greater than a tiny threshold
            // (this avoids division by zero for silent bands)
            if (maxEnergy > 1e-9f) {
                for (float& energy_value : band_envelope) {
                    energy_value /= maxEnergy;
                }
            }
        }
    }
}