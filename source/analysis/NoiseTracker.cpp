#include "NoiseTracker.h"

NoiseTracker::NoiseTracker(const Sitrano::NoiseSettings& settings,
    const Sitrano::AnalysisUnit& ana, 
	const Sitrano::Results& r,
	float startFreq):
	unit(ana),
    sf(unit.soundFile),
    N(settings.nfft),
    topFreqs(r.topFreqs),
    hop(settings.hopSize),
    num(settings.numBins),
    sr(unit.sampleRate),
    useList(settings.useList),
	startFreq(startFreq)
{
	float f0 = settings.minFreq;
    float m0 = Sitrano::freqToMidi(f0);
    float step;
    float half;
    if (settings.useOctaveDiv) step = 12.0 / settings.octaveDiv;
    else step = (Sitrano::freqToMidi(settings.maxFreq) - Sitrano::freqToMidi(f0)) / 
        (float)num;
    half = step / 2.0;
	while (m0 < Sitrano::freqToMidi(sr / 2.0) && 
        m0 < Sitrano::freqToMidi(settings.maxFreq)) 
    {
		float f_low = Sitrano::midiToFreq(m0 - half);
        float f_high = Sitrano::midiToFreq(m0 + half);
		bands.push_back({ f_low, f_high, f0 });
        m0 += step;
		f0 = Sitrano::midiToFreq(m0);
	}

	// FFTW Initialization
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
    useList = true;
}

std::vector<std::vector<float>> NoiseTracker::analyze() {
    // --- Initialization (mostly unchanged) ---
    const std::vector<float>& signal = unit.soundFile;
    const float sampleRate = unit.sampleRate;
    std::vector<std::vector<float>> noise;

    int numFrames = useList ? frameTable.size() : (signal.size() - N) / hop;
    if (numFrames <= 0) return noise;
    const float df = sampleRate / N;

    // The outer vector now represents BANDS, the inner vector represents FRAMES (time).
    noise.assign(bands.size(), std::vector<float>(numFrames, 0.0f));

    // Main Analysis Loop (calculates power)
    for (int i = 0; i < numFrames; ++i) {
        int startSample = useList ? frameTable[i] : i * hop;
        if (startSample + N > signal.size()) {
            // Trim all band vectors if the signal ends unexpectedly
            for (auto& band_vec : noise) {
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
            float mag = real * real + imag * imag;
            float amp = Sitrano::mag_to_amp(mag, N);
            int bin = findLargeBin(freq, bands);

            noise[bin][i] += amp;
        }
    }

    // --- 2. Normalization Step ---
    // This runs after all power values have been calculated for all frames.
    //for (auto& band_envelope : noise) {
    //    // Find the peak energy in this band's envelope using std::max_element
    //    auto max_iterator = std::max_element(band_envelope.begin(), band_envelope.end());

    //    // Proceed only if the band is not empty
    //    if (max_iterator != band_envelope.end()) {
    //        float maxEnergy = *max_iterator;

    //        // Normalize the envelope if the peak is greater than a tiny threshold
    //        // (this avoids division by zero for silent bands)
    //        if (maxEnergy > 1e-9f) {
    //            for (float& energy_value : band_envelope) {
    //                energy_value /= maxEnergy;
    //            }
    //        }
    //    }
    //}
    return noise;
}

bool NoiseTracker::isInVector(float freq, const std::vector<Sitrano::Peak>& topFreqs)
{
    bool isInVector = false;
    int freqBin = Sitrano::freqToBin(freq, N, sr);
    for (int i = 0; i < topFreqs.size(); i++)
    {
        int topFreqBin = Sitrano::freqToBin(topFreqs[i].freq, N, sr);
        if (topFreqBin == freqBin) isInVector = true;
    }

    return isInVector;
}

int NoiseTracker::findLargeBin(float freq, std::vector<Band> bands)
{
    int largeBin = 0;
    for (int i = 0; i < bands.size(); i++)
    {
        if (freq > bands[i].f_low && freq <= bands[i].f_high)
        {
            largeBin = i;
            return largeBin;
        }
    }
    return largeBin;
}