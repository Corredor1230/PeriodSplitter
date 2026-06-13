#pragma once

#include<vector>
#include<string>
#include<cmath>
#include<cstdint>
#include<type_traits>
#include<numeric>
#include<iostream>
#include<algorithm>
#include<complex>

#include"ResynthHeader.h"

namespace Sihat
{
    constexpr double PI = 3.14159265358979323846;
    constexpr double M_E = 2.718281828459045;
    constexpr double TWO_PI = 2.0 * PI;
	//Structures
    // Enum for different FFT window styles used for analysis
    enum WindowStyle {
        periodLoop = 0,
        singlePeriod = 1,
        audioChunk = 2
    };

    enum class WaveletType {
        morlet = 0,
        meyer = 1,
        haar = 2,
        hat = 3
    };
    
    struct HarmonicResults {
        std::vector<std::vector<float>> amps;
        std::vector<std::vector<float>> phases;
        std::vector<std::vector<float>> freqs;
        std::vector<uint32_t> finalSamples;
        float rms = 0.0;
    };

    struct ComplexSpectrogram {
        std::vector<std::complex<float>> data;
        int numFrames;
        int numBins;
        // A helper to make access look clean again
        std::complex<float>& at(int frame, int bin) {
            return data[frame * numBins + bin];
        }

        const std::complex<float>& at(int frame, int bin) const {
            return data[frame * numBins + bin];
        }
    };

    struct Spectrogram{
        std::vector<float> data;
        int numFrames;
        int numBins;
        float& at(int frame, int bin){
            return data[frame * numBins + bin];
        }

        const float& at(int frame, int bin) const {
            return data[frame * numBins + bin];
        }
    };

    struct VariableRatePartial
    {
        float frequency;
        int hopSize;
        std::vector<float> data; // The envelope, decimated
    };

    // Main settings for the Analyzer. 
    // This structure decides which processing steps will be completed.
    struct Settings {
        bool fullHarmonicAnalysis = true;
        bool sourceSeparation = true;
        bool pitchAnalysis = true;
        bool transientSeparation = true;
        bool periodAnalysis = true;
        bool overtoneAnalysis = true;
        bool transientAnalysis = true;
        bool harmonicAnalysis = true;
        bool noiseAnalysis = true;
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
        int tailOff = 500;
        float tCorrelationThreshold = 0.85;
    };

    struct TransientFFTSettings {
        bool useFFT = true;
        int nfft = 1024;
        int hopSize = 256;
        float flatnessThreshold = 0.6;
    };

    struct TWaveletSettings {
        WaveletType wtype = WaveletType::morlet;
        bool forceFFTSize = false;
        int fftSize = 4096;
    };

    struct SingleTransientSettings{
        int nfft = 512;
        int hopSize = 32;
        int rmsWindow = 30;
        int rmsHopSize = 15;
        int startSample = 0;
        int numBands = 16;
        int maxOvertones = 32;
        int overFrames = 16;
        int overNfft = 2048;
        float outThreshold = 0.05f;
        float inThreshold = 0.1f;
        float tolInCents = 50.0;
        float o_tolInCents = 10.0;
        int numModes = 32;
        int modeWindow = 2048;
        float firstPeakThreshold = 0.4f;
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
        bool applyHanning = true;
        WindowStyle style = WindowStyle::audioChunk;
        float toleranceValue = 100.0;
        bool useTolInHz = true;
        float tolInHz = 6.0;
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
        bool useTolInHz = true;
        double tolInHz = 6.0;
        bool postTransientStart = true;
        bool chooseFirstSample = false;
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

    struct ExportSettings {
        bool exportHarmonicPhase = false;
        bool exportSpectralCentroid = false;
        bool exportSpectralFlatness = true;
    };

    struct NoiseSettings {
        int nfft = 128;
        int hopSize = 32;
        int numBins = 32;
        float minFreq = 60.f;
        float maxFreq = 20000.f;
        bool useOctaveDiv = true;
        float octaveDiv = 3.0;
        int startSample = 0;
        bool useList = false;
    };

    struct HPSSSettings {
        int nfft = 4096;
        int hopSize = 1024;
        int filtSize = 17;
        float exponent = 2.0;
    };

    struct AnalysisConfig {
        int numHarmonics = 32;
        int nfft = 16384 * 2;
        int hopSize = 1024;
        int startSample = 0;
        float tolerance = 100.f;
        std::string outDir = "../../media/tests";
        PitchSettings pSettings;
        HPSSSettings hpSettings;
        TransientSettings tSettings;
        TransientFFTSettings tfftSettings;
        SingleTransientSettings stSettings;
        CorrelationSettings cSettings;
        OvertoneSettings oSettings;
        HarmonicSettings hSettings;
        NoiseSettings nSettings;
        bool verbose = true;
        bool bulkProcess = false;
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
        double pha = 0.0;
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
        Sihat::BinFreq bin;
        float mag;
        float amp;
        float pha;
    };

    struct TrackedPoint {
        float freq;          // The exact interpolated frequency at this frame
        float crestFactor;   // Prominence vs noise floor (can be swapped for absolute amp/mag if preferred later)
        float amp;
        bool active;          // True if a distinct peak was found within tolerance
    };

    //The information for a single transient overtone. It contains the target and a vector with frames.
    struct TransientOvertone {
        Sihat::Peak targetOvertone; 
        std::vector<TrackedPoint> envelope;
    };

    //Outdated, eventually will be fully deprecated.
    struct TransientResults {
        SampleRange range;
        std::vector<VariableRatePartial> scalogram;
        float rms = 0.0;
    };

    /**
     * @brief This structure contains the more resonant aspects of a transient. Therefore, it is modeled similarly to the Harmonic Analysis. 
     * @param overtones: A vector of overtones. Each of them will contain a target overtone as well as an envelope.
     * @param floor: A vector of floor values. These are simple amplitude averages of each frame's spectrum. Useful to calculate SNR values for overtones or other peaks.
    */
    struct TransientHarmonics{
        uint32_t hopSize = 16;
        uint32_t startSample = 0;
        //Each overtone contains one target and one envelope.
        std::vector<TransientOvertone> overtones;
        //Should have one floor per frame. So the size should be equal to overtones[n].envelope.size()
        std::vector<float> floor;
    };

    struct STransientResults{
        SampleRange range;
        uint32_t envHopSize;
        uint32_t specHopSize;
        uint32_t floorHopSize;
        uint32_t specWindowSize;
        uint32_t specNumBins;
        uint32_t specFrameNum;
        int riseTime;
        float peakAmp;
        std::vector<Sample> ampEnvelope;
        std::vector<float> centroid;
        std::vector<float> flatness;
        std::vector<float> bandEnvelopes;
        Synth::TModes tmodes;
        TransientHarmonics tHarmonics;
        float rms;
    };



	//Larger structure containing all results of the analysis
	struct Results {
        TransientResults tResults;
        STransientResults stResults;
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
    inline float minmax(float low, float high, float in)
    {
        return std::max<float>(std::min<float>(high, in), low);
    }
    inline std::vector<float> getAudioExcerpt(const std::vector<float>& audio, int start, int end)
    {
        if (start < 0 || end > audio.size() || start >= end) {
            return {};
        }
        return std::vector<float>(audio.begin() + start, audio.begin() + end);
    }
    inline std::vector<float> getDifferential(const std::vector<float>& audio)
    {
        std::vector<float> diff(audio.size(), 0.0f);
        for (size_t i = 1; i < audio.size(); ++i)
        {
            diff[i] = audio[i] - audio[i - 1];
        }
        return diff;
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
    inline int findNextPowerOfTwo(int size)
    {
        if (size < 0)
        {
            std::cerr << "Error: Non-positive size" << '\n';
        }
        for (int i = 1; i < 32; i++)
        {
            int prevPow = static_cast<int>(std::pow(2.0, static_cast<float>(i) - 1));
            int currPow = static_cast<int>(std::pow(2.0, static_cast<float>(i)));

            if (size > prevPow && size < currPow) return currPow;
        }

        std::cerr << "Error: Transient size too big" << '\n';
        return 0;
    }
    inline float findMaxAbsVal(const std::vector<float>& in)
    {
        float max = 0.0;
        for (int i = 0; i < in.size(); i++)
        {
            if (std::fabs(in[i]) > max)
            {
                max = in[i];
            }
        }
        return max;
    }
    inline float normLogistic(float in, float exp)
    {
        float lin = std::clamp(in, 0.0f, 1.0f);
    
        // 2. Exact boundary handling to prevent NaN and 0/0 errors
        if (lin <= 0.0f) return 0.0f;
        if (lin >= 1.0f) return 1.0f;
        
        // 3. Handle exponent edge cases (fallback to linear or step)
        if (exp <= 0.0f) return 0.5f; // Flat mid-point if exponent is invalid/zero
        if (exp == 1.0f) return lin;  // Perfect linear shortcut
        
        // 4. Core calculation (now safely away from 0.0 and 1.0)
        float num = std::pow(lin, exp);
        float den = num + std::pow(1.0f - lin, exp);
        
        return num / den;
    }
    inline int findPrevPowerOfTwo(int size)
    {
        if (size < 0)
        {
            std::cerr << "Error: non-positive size" << '\n';
        }
        for (int i = 1; i < 32; i++)
        {
            int prevPow = static_cast<int>(std::pow(2.0, static_cast<float>(i) - 1));
            int currPow = static_cast<int>(std::pow(2.0, static_cast<float>(i)));

            if (size > prevPow && size < currPow) return currPow;
        }

        std::cerr << "Error: size too big" << '\n';
        return 0;
    }
    inline float getRmsValue(const std::vector<float>& input, int start, int end)
    {
        if (start >= end || end > input.size() || input.empty()) {
            return 0.0;
        }

        double sumOfSquares = 0.0;
        std::for_each(input.begin() + start, input.begin() + end, [&sumOfSquares](double val) {
            sumOfSquares += val * val;
        });

        size_t count = end - start;
        
        // Return RMS: sqrt(sum of squares / number of elements)
        return std::sqrt(sumOfSquares / count);
    }
    inline float getMaxVal(const std::vector<float>& input)
    {
        float maxval = 0.0;
        for (int i = 0; i < input.size(); i++)
        {
            if (input[i] > maxval) maxval = input[i];
        }

        return maxval;
    }

    //Functions with a definition in .cpp file
    BinFreq findPeakWithinTolerance(float targetFreq, float tolerance, int n, float sr, void* out, bool isTolInHz);
    FreqUnit findPeak(BinFreq target, void* fftwfOut, int nfft, float fs);
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
        const HarmonicResults& hResults,
        const TransientResults& tResults,
        float f0,
        const std::string& outputDirectory,
        const std::string& baseFilename,
        const std::string& extension
    );
    void saveHarmonicDataSihat(
        const HarmonicResults& hResults,
        const STransientResults& stResults,
        const AnalysisConfig& config,
        const Settings& settings,
        const float f0,
        const uint32_t sr,
        const std::string& outputDirectory,
        const std::string& prefix,
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

    inline std::vector<float> getHannWindow(int size)
    {
        std::vector<float> window(size);
        for(int i = 0; i < size; i++)
        {
            window[i] = 0.5f * (1.0f - std::cos(2.0f * PI * i / (size - 1.0f)));
        }

        return window;
    }

    inline float getHannValue(int i, int nfft)
    {
        if (i < 0 || i >= nfft) return 0.0;
        float windowVal = 0.5 * (1.0 - std::cos(TWO_PI * i / (double)(nfft - 1)));
        return windowVal;
    }



};