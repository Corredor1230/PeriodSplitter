#pragma once

#include<iostream>
#include<vector>
#include<string>
#define _USE_MATH_DEFINES
#include<math.h>
#include<cmath>
#include<fstream>
#include<cstdint>
#include<fftw3.h>
#include<filesystem>
#include<type_traits>
#include<iomanip>

namespace Sitrano
{

    constexpr double PI = 3.14159265358979323846;
	//Structures
    // Enum for different FFT window styles used for analysis
    enum WindowStyle {
        periodLoop = 0,
        singlePeriod,
        audioChunk
    };

    struct HarmonicResults {
        std::vector<std::vector<float>> amps;
        std::vector<std::vector<float>> phases;
        std::vector<std::vector<float>> freqs;
        std::vector<uint32_t> finalSamples;
    };

    /* This structure contains the main settings for the HarmonicTracker class
    * applyHanning: Defines whether a hanning window is applied to each FFT window.
    * windowStyle: Features three different types of FFT processing
    |-----periodLoop: Fills the FFT window with a single period looped over.
    |-----singlePeriod: zero-pads the window after the single harmonci event.
    L-----audioChunk: takes the window directly from the original audio signal.
    * toleranceValue: The value in cents around each top harmonic where buckets will be considered the same overtone.
    */
    struct HarmonicSettings {
        bool applyHanning;
        WindowStyle style;
        float toleranceValue;
    };

    struct OvertoneSettings {
        bool useTolerance = true;
        int overtoneFirstSample = 0;
        double toleranceValue = 100.0;
        bool useCustomSignal = true;
        bool sumAmplitudes = true;
        int fftSize = 65536;
    };

    struct AnalysisConfig {
        int numHarmonics = 32;
        int nfft = 16384 * 2;
        int hopSize = 1024;
        int startSample = 0;
        float tolerance = 100.f;
        OvertoneSettings oConfig{
            true, 0, 100.0, true, true, 65536
        };
        HarmonicSettings hConfig{
            true, WindowStyle::audioChunk, 100.f
        };
    };

    struct ChangePoint {
        uint32_t index; // 4 bytes
        float value;    // 4 bytes
        float ratio;
    };

    struct FileHeader {
        uint32_t numChannels;
        uint32_t numSamples;
        float ratioConstant;
    };

	//This contains the necessary info for all analysis classes
	struct AnalysisUnit {
		std::vector<float> soundFile;
        std::string filename;
		float sampleRate;
	};

	//An amplitude bucket resulting from an FFT
	struct Peak {
		double freq = 0.0;
		double mag = 0.0;
		double amp = 0.0;
	};

    // Main settings for the Analyzer. 
    // This structure decides which processing steps will be completed.
    struct Settings {
        bool pitchAnalysis = true;
        bool periodAnalysis = true;
        bool overtoneAnalysis = true;
        bool harmonicAnalysis = true;
        bool noiseAnalysis = true;
        bool transientSeparation = true;
    };

	//Larger structure containing all results of the analysis
	struct Results {
        std::vector<uint32_t> sampleList;
		std::vector<Peak> topFreqs;
        float pitch;
        std::vector<std::vector<float>> noise;
        HarmonicResults hResults;
	};

	//Helper functions
    static inline double hann_gain_rms() { return 0.5; };
    static inline double interp_delta(int k, const std::vector<double>& mags)
    {
        if (k <= 0 || k >= (int)mags.size() - 1) return 0.0;
        double m1 = mags[k - 1], m0 = mags[k], p1 = mags[k + 1];
        double denom = (m1 - 2 * m0 + p1);
        if (std::abs(denom) < 1e-20) return 0.0;
        return 0.5 * (m1 - p1) / denom;
    }
    static inline double mag_to_amp(double mag, int N)
    {
        double A = (2.0 * mag) / double(N);
        A /= hann_gain_rms();
        return A;
    } //Converts FFT magnitude to linear gain
    static inline float dbToAmp(float db)
    {
        float gain = std::pow<float>(10.0, db / 20.0);
        return gain;
    } //Converts FS decibel values to linear gain
    static inline float ampToDb(float amp)
    {
        float decibels = 20.0 * std::log10(amp / 1.0);
        return decibels;
    } //Converts linear gain to FS decibel
    static inline float midiToFreq(float midi) {
        float freq = std::pow<float>(2.0, (midi - 69.0) / 12.0) * 440.0;
        return freq;
    }
    static inline float freqToMidi(float freq) {
        float midi = 12.0 * std::log2(freq / 440.0) + 69.0;
        return midi;
    }
    static inline double cents_to_hz(double f_center, double cents)
    {
        double r = std::pow<double>(2.0, cents / 1200.0);
        return std::max<double>(1e-12, f_center * (r - 1.0));
    }
    static inline int freqToBin(float freq, int n, float fs) {
        int bin = static_cast<int>(std::round(freq * n / fs));
        return bin;
    }
    static inline float binToFreq(int bin, int n, float fs) {
        float freq = n != 0 ? static_cast<float>(bin) * fs / static_cast<float>(n) : 0.0;
        return freq;
    }
    static inline float findPeakWithinTolerance(float targetFreq, float tolerance, int n, float sr, fftwf_complex* out) {
        float midiTolerance = tolerance / 100.0;
        float lowMidi = freqToMidi(targetFreq) - midiTolerance;
        float lowFreq = midiToFreq(lowMidi);
        int lowBin = freqToBin(lowFreq, n, sr);
        float hiMidi = freqToMidi(targetFreq) + midiTolerance;
        float hiFreq = midiToFreq(hiMidi);
        int hiBin = freqToBin(hiFreq, n, sr);
        int targetBin = freqToBin(targetFreq, n, sr);


        int binNumber = std::abs(hiBin - lowBin);

        float outFreq = 0.f;
        float peakAmp = 0.f;
        if (binNumber != 0)
        {
            for (int i = 0; i < binNumber; i++)
            {
                int currentBin = lowBin + i;
                float real = out[currentBin][0];
                float imag= out[currentBin][0];

                float mag = std::sqrt(imag * imag + real * real);
                float amp = mag_to_amp(mag, n);

                if (amp > peakAmp)
                {
                    outFreq = binToFreq(currentBin, n, sr);
                    peakAmp = amp;
                }
                
            }
        }
        else outFreq = targetFreq;

        return outFreq;
    };
    static inline void filterVector(std::vector<float>& input, int filterSize, bool zeroPad) {
        std::vector<float> filter;
        filter.resize(filterSize);
        if (zeroPad)
        {
            for (int i = 0; i < filter.size(); i++)
            {
                filter[i] = 0.f;
            }
        }
        else
        {
            for (int i = 0; i < filter.size(); i++)
            {
                filter[i] = input[0];
            }
        }

        for (int i = 0; i < input.size(); i++)
        {
            int index = i % filterSize;
            filter[index] = input[i];
            float sum = 0.f;
            for (int j = 0; j < filter.size(); j++)
            {
                sum += filter[j];
            }
            if (filter.size() != 0) sum /= (float)filter.size();
            input[i] = sum;
        }
    }
    static inline void normalizeByMaxAbs(std::vector<float>& vec)
    {
        if (vec.empty()) return;

        // Find the maximum absolute value
        float maxAbs = 0.0f;
        for (float x : vec) {
            maxAbs = std::max<float>(maxAbs, std::abs(x));
        }

        if (maxAbs == 0.0f) {
            // Avoid division by zero — all elements are zero
            return;
        }

        // Normalize
        for (float& x : vec) {
            x /= maxAbs;
        }
    }
    static inline std::string getRawFilename(std::string& filename)
    {
        std::string rawFilename;
        std::filesystem::path path = filename;
        rawFilename = path.filename().replace_extension("").string();

        return rawFilename;
    }
    static inline void applyWindow(std::vector<float>& frame) {
        int size = frame.size();
        for (int i = 0; i < size; ++i) {
            frame[i] *= 0.5f * (1.0f - cosf(2.0f * M_PI * i / (size - 1)));
        }
    }

    static float getPitchFromFilename(const std::string& filename)
    {
        std::string name = filename;
        std::vector<std::string> parts;
        size_t start = 0, pos;

        while ((pos = name.find('_', start)) != std::string::npos) {
            parts.push_back(name.substr(start, pos - start));
            start = pos + 1;
        }
        parts.push_back(name.substr(start));

        float outPitch = std::stof(parts[2]);

        return outPitch;
    }

    static void saveHarmonicData(
        const std::vector<uint32_t>& indices,
        const std::vector<std::vector<float>>& fastData,   // RMS per harmonic
        const std::vector<std::vector<float>>& slowData,   // raw frequencies per harmonic
        float fundamentalFreq,                             // reference frequency
        const std::string& pathIndex,
        const std::string& pathFast,
        const std::string& pathSlow
    ) {
        // Write indices
        {
            std::ofstream f(pathIndex, std::ios::binary);
            uint32_t n = static_cast<uint32_t>(indices.size());
            f.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
            f.write(reinterpret_cast<const char*>(indices.data()), n * sizeof(uint32_t));
        }

        // Write fast RMS data
        {
            std::ofstream f(pathFast, std::ios::binary);
            uint32_t numHarmonics = static_cast<uint32_t>(fastData.size());
            f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));
            for (const auto& harmonic : fastData) {
                uint32_t len = static_cast<uint32_t>(harmonic.size());
                f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(harmonic.data()), len * sizeof(float));
            }
        }

        // Write slow frequency data (auto-generate change points)
        {
            std::ofstream f(pathSlow, std::ios::binary);
            uint32_t numHarmonics = static_cast<uint32_t>(slowData.size());
            f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));

            for (const auto& freqVec : slowData) {
                std::vector<ChangePoint> cps;

                if (!freqVec.empty()) {
                    float prev = freqVec.front();
                    cps.push_back({ indices.front(), prev, prev / fundamentalFreq });

                    for (size_t i = 1; i < freqVec.size() && i < indices.size(); ++i) {
                        if (std::fabs(freqVec[i] - prev) > 1e-6f) { // detect change
                            cps.push_back({ indices[i], freqVec[i], freqVec[i] / fundamentalFreq });
                            prev = freqVec[i];
                        }
                    }
                }

                uint32_t len = static_cast<uint32_t>(cps.size());
                f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                for (const auto& cp : cps) {
                    f.write(reinterpret_cast<const char*>(&cp.index), sizeof(uint32_t));
                    f.write(reinterpret_cast<const char*>(&cp.value), sizeof(float));
                    f.write(reinterpret_cast<const char*>(&cp.ratio), sizeof(float));
                }
            }
        }
    }

    template <typename T>
    static void applyHann(T* data, size_t N) {
        static_assert(std::is_arithmetic_v<T>,
            "applyHann can only be used with numerical types.");

        // 1. Define constants for precision
        constexpr double TWO_PI = 2.0 * M_PI;

        // 2. Handle edge cases
        if (N < 2) {
            if (N == 1) {
                data[0] = static_cast<T>(0.0);
            }
            return; // Do nothing for size 0
        }

        // 3. Denominator for window (N-1)
        // We use 'double' for all intermediate math to maintain precision.
        const double M = static_cast<double>(N - 1);

        // 4. Apply the window
        for (size_t n = 0; n < N; ++n) {
            const double n_double = static_cast<double>(n);

            // Hanning window formula: w(n) = 0.5 * (1 - cos(2*pi*n / (N-1)))
            const double window_val = 0.5 * (1.0 - std::cos(TWO_PI * n_double / M));

            // 5. Apply the window in-place, casting the final result back to T
            data[n] = static_cast<T>(static_cast<double>(data[n]) * window_val);
        }
    }

    template <typename T>
    static void applyHann(std::vector<double>& data) {
        static_assert(std::is_arithmetic_v<T>,
            "applyHanningWindow can only be used with numerical types.");

        const size_t N = data.size();

        // 2. Handle edge cases to avoid division by zero
        if (N < 2) {
            if (N == 1) {
                data[0] = static_cast<T>(0.0); // A 1-point Hanning window is 0
            }
            return; // Do nothing for an empty vector
        }

        // 3. We use 'double' for all intermediate math to maintain precision,
        //    regardless of whether T is float, double, or int.
        const double M = static_cast<double>(N - 1);
        constexpr double TWO_PI = 2.0 * PI;

        for (size_t n = 0; n < N; ++n) {
            const double n_double = static_cast<double>(n);

            // Hanning window formula: w(n) = 0.5 * (1 - cos(2*pi*n / (N-1)))
            const double window_val = 0.5 * (1.0 - std::cos(TWO_PI * n_double / M));

            // 4. Apply the window. We explicitly cast the final result
            data[n] = static_cast<T>(data[n] * window_val);
        }
    }



};