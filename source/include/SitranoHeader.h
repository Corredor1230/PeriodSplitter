#pragma once

#include<vector>
#include<string>
#include<cmath>
#include<cstdint>
#include<type_traits>

namespace Sitrano
{

    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
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

    struct TransientSettings {
        int tStartSample = 0;
        bool useMs = false;
        int rmsSampleSize = 128;
        int rmsSampleHopLength = 64;
        float transientRmsSizeMs = 1.0f;
        float transientRmsHopRatio = 1.0f;
        float transientFactor = 3.0f;
        float transientThreshold = 0.1f;
        int preAttack = 20; // number of samples before transient detection
    };

    /** 
    * @brief This structure contains the main settings for the HarmonicTracker class
    * @param applyHanning: Defines whether a hanning window is applied to each FFT window.
    * @param windowStyle: Features three different types of FFT processing
    * @param periodLoop: Fills the FFT window with a single period looped over.
    * @param singlePeriod: zero-pads the window after the single harmonci event.
    * @param audioChunk: takes the window directly from the original audio signal.
    * @param toleranceValue: The value in cents around each top harmonic where buckets will be considered the same overtone.
    */
    struct HarmonicSettings {
        bool applyHanning;
        WindowStyle style;
        float toleranceValue;
    };

    /**
    * @brief This structure contains the main parameters for the most representative overtone analysis.
    * @param useTolerance: Decides whether bins should be grouped as one.
    * @param toleranceValue: Measured in cents
    * @param overtoneFirstSample: The first sample from which the top overtones will be found.
    * @param useCustomSignal: If false, no custom data range will be analyzed.
    * @param sumAmplitudes: Sums the amplitudes of all bins within the tolerance range.
    * @param fftSize: Generally a larger FFT size than the global config.
    * toleranceValue is measured in cents.
    */
    struct OvertoneSettings {
        bool useTolerance = true;
        double toleranceValue = 100.0;
        bool postTransientStart = true;
        int overtoneFirstSample = 0;
        bool useCustomSignal = true;
        bool sumAmplitudes = true;
        int fftSize = 65536;
        float overtoneThreshold = -60.f;
        bool setAbsoluteThreshold = true;
    };
    struct PitchSettings {
        float modeThreshold = 3.0;
        float toleranceInCents = 50.f;
        float minFreq = 60.f;
        float maxFreq = 1300.f;
    };

    struct CorrelationSettings {
        //Period detection
        float periodStartOffsetMs = 50.0f;
        float correlationThreshold = 0.95f;
    };

    struct AnalysisConfig {
        int numHarmonics = 32;
        int nfft = 16384 * 2;
        int hopSize = 1024;
        int startSample = 0;
        float tolerance = 100.f;
        PitchSettings pConfig;
        TransientSettings tSettings;
        CorrelationSettings cSettings;
        OvertoneSettings oConfig;
        HarmonicSettings hConfig;
        bool verbose = true;
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

    struct Sample {
        int index;
        float value;
    };

    struct BinFreq {
        int bin;
        float freq;
        float fBin;
    };

    struct SampleRange {
        int initSample;
        int endSample;
    };

    struct FreqUnit {
        Sitrano::BinFreq bin;
        float mag;
        float amp;
        float pha;
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
        SampleRange tRange;
        std::vector<uint32_t> sampleList;
		std::vector<Peak> topFreqs;
        float pitch;
        std::vector<std::vector<float>> noise;
        HarmonicResults hResults;
	};

	//Small helper functions
    inline double hann_gain_rms() { return 0.5; };
    inline double interp_delta(int k, const std::vector<double>& mags)
    {
        if (k <= 0 || k >= (int)mags.size() - 1) return 0.0;
        double m1 = mags[k - 1], m0 = mags[k], p1 = mags[k + 1];
        double denom = (m1 - 2 * m0 + p1);
        if (std::abs(denom) < 1e-20) return 0.0;
        return 0.5 * (m1 - p1) / denom;
    }
    inline double mag_to_amp(double mag, int N)
    {
        double A = (2.0 * mag) / double(N);
        A /= hann_gain_rms();
        return A;
    } //Converts FFT magnitude to linear gain
    /**
    * @brief Converts linear amplitude to FFT magnitude for a non-DC/non-Nyquist bin.
    *
    * This is the inverse of the mag_to_amp function.
    *
    * @param A The linear amplitude of the sinusoid.
    * @param N The FFT size.
    * @return double The corresponding FFT magnitude bin value.
    */
    inline double amp_to_mag(double A, int N)
    {
        double mag = (A * (double)N) / 4.0;
        return mag;
    }
    inline float dbToAmp(float db)
    {
        float gain = std::pow<float>(10.0, db / 20.0);
        return gain;
    } //Converts FS decibel values to linear gain
    inline float ampToDb(float amp)
    {
        float decibels = 20.0 * std::log10(amp / 1.0);
        return decibels;
    } //Converts linear gain to FS decibel
    inline float midiToFreq(float midi) {
        float freq = std::pow<float>(2.0, (midi - 69.0) / 12.0) * 440.0;
        return freq;
    }
    inline float freqToMidi(float freq) {
        float midi = 12.0 * std::log2(freq / 440.0) + 69.0;
        return midi;
    }
    inline double cents_to_hz(double f_center, double cents)
    {
        double r = std::pow<double>(2.0, cents / 1200.0);
        return std::max<double>(1e-12, f_center * (r - 1.0));
    }
    inline int freqToBin(float freq, int n, float fs) {
        int bin = static_cast<int>(std::round(freq * n / fs));
        return bin;
    }
    inline float binToFreq(int bin, int n, float fs) {
        float freq = n != 0 ? static_cast<float>(bin) * fs / static_cast<float>(n) : 0.0;
        return freq;
    }
    inline float binToFreq(float bin, int n, float fs)
    {
        float freq = n != 0 ? bin * fs / static_cast<float>(n) : 0.0;
        return freq;
    }
    inline bool withinTolerance(float source, float comp, float tolerance)
    {
        bool result = false;
        if (std::abs(source - comp) < tolerance) result = true;
        return result;
    }

    //Functions with a definition in .cpp file
    BinFreq findPeakWithinTolerance(float targetFreq, float tolerance, int n, float sr, void* out);
    FreqUnit findPeak(BinFreq target, void* fftwfOut, int nfft, float fs, int binRange);
    void filterVector(std::vector<float>& input, int filterSize, bool zeroPad);
    void normalizeByMaxAbs(std::vector<float>& vec);
    std::string getRawFilename(const std::string& filename);
    void applyWindow(std::vector<float>& frame);
    float getPitchFromFilename(const std::string& filename);
    void saveHarmonicData(
        const std::vector<uint32_t>& indices,
        const std::vector<std::vector<float>>& fastData,   // RMS per harmonic
        const std::vector<std::vector<float>>& slowData,   // raw frequencies per harmonic
        float fundamentalFreq,                             // reference frequency
        const std::string& pathIndex,
        const std::string& pathFast,
        const std::string& pathSlow
    );
    void saveHarmonicDataSihat(
        const std::vector<uint32_t>& indices,
        const std::vector<std::vector<float>>& amp,
        const std::vector<std::vector<float>>& freq,
        float f0,
        const std::string& outputDirectory,
        const std::string& baseFilename,
        const std::string& extension
    );
    int findPeakIndexVector(const std::vector<float>& input);
    int findPreviousZero(const std::vector<float>& signal, int startSmaple);
    int findNextZero(const std::vector<float>& signal, int startSample);
    int findNearestZero(const std::vector<float>& signal, int startSample);
    int findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute = true);
    std::vector<int> findZeroCrossings(const std::vector<float>& signal, int startSample);
    int findNearestCachedZero(const std::vector<int>& zeroCrossings, int sample);
    int findAbsPeakIndex(const std::vector<float>& vector);


    //Template functions
    template <typename T>
    static void applyHann(T* data, size_t N) {
        static_assert(std::is_arithmetic_v<T>,
            "applyHann can only be used with numerical types.");

        // Handle edge cases
        if (N < 2) {
            if (N == 1) {
                data[0] = static_cast<T>(0.0);
            }
            return;
        }

        // We use 'double' for all intermediate math to maintain precision.
        const double M = static_cast<double>(N - 1);

        // Apply the window
        for (size_t n = 0; n < N; ++n) {
            const double n_double = static_cast<double>(n);
            const double window_val = 0.5 * (1.0 - std::cos(TWO_PI * n_double / M));

            // Apply the window in-place, casting the final result back to T
            data[n] = static_cast<T>(static_cast<double>(data[n]) * window_val);
        }
    }

    template <typename T>
    static void applyHann(std::vector<double>& data) {
        static_assert(std::is_arithmetic_v<T>,
            "applyHanningWindow can only be used with numerical types.");

        const size_t N = data.size();

        // Handle edge cases to avoid division by zero
        if (N < 2) {
            if (N == 1) {
                data[0] = static_cast<T>(0.0); // A 1-point Hanning window is 0
            }
            return;
        }

        // We use double for all intermediate math to maintain precision
        const double M = static_cast<double>(N - 1);
        constexpr double TWO_PI = 2.0 * PI;

        for (size_t n = 0; n < N; ++n) {
            const double n_double = static_cast<double>(n);
            const double window_val = 0.5 * (1.0 - std::cos(TWO_PI * n_double / M));

            // We explicitly cast the final result
            data[n] = static_cast<T>(data[n] * window_val);
        }
    }



};