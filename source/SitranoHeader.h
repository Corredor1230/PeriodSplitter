#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<cmath>
#include<fftw3.h>

class Sitrano
{
public:
	//Structures

	//This contains the necessary info for all analysis classes
	struct AnalysisUnit {
		std::vector<float> soundFile;
		std::vector<int> sampleList;
		float sampleRate;
		float pitch;
		int numHarmonics;
		int nfft;
	};

	//An amplitude bucket resulting from an FFT
	struct Peak {
		double freq = 0.0;
		double mag = 0.0;
		double amp = 0.0;
	};

	//Larger structure containing all results of the analysis
	struct Results {
		std::vector<std::vector<float>> amps;
		std::vector<std::vector<float>> phases;
		std::vector<std::vector<float>> freqs;
		std::vector<Peak> topFreqs;
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
        float midi = 12.0 * std::log2(freq / 440.0);
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
        if (binNumber != 0)
        {
            for (int i = 0; i < binNumber - 1; i++)
            {
                int currentBin = lowBin + i;
                float nextReal = out[currentBin + 1][0];
                float prevReal = out[currentBin][0];
                float nextImag = out[currentBin + 1][1];
                float prevImag = out[currentBin][0];

                float nextMag = std::sqrt(nextImag * nextImag + nextReal * nextReal);
                float nextAmp = mag_to_amp(nextMag, n);

                float prevMag = std::sqrt(prevImag * prevImag + prevReal * prevReal);
                float prevAmp = mag_to_amp(prevMag, n);

                prevAmp < nextMag ?
                    outFreq = binToFreq(currentBin, n, sr) :
                    outFreq = binToFreq(currentBin + 1, n, sr);

            }
        }
        else outFreq = targetFreq;

        return outFreq;
    };
};