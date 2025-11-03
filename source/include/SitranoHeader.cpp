#include "SitranoHeader.h"
#include<iostream>
#include<fftw3.h>
#include<fstream>
#include<filesystem>
#include<iomanip>
#include<complex>

Sitrano::BinFreq Sitrano::findPeakWithinTolerance(float targetFreq, float tolerance, int n, float sr, void* generic_out)
{
    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(generic_out);
    float midiTolerance = tolerance / 100.0;
    float lowMidi = freqToMidi(targetFreq) - midiTolerance;
    float lowFreq = midiToFreq(lowMidi);
    int lowBin = freqToBin(lowFreq, n, sr);
    float hiMidi = freqToMidi(targetFreq) + midiTolerance;
    float hiFreq = midiToFreq(hiMidi);
    int hiBin = freqToBin(hiFreq, n, sr);
    int targetBin = freqToBin(targetFreq, n, sr);

    Sitrano::BinFreq outBF{ 0.f, 0.f };

    int binNumber = std::abs(hiBin - lowBin);

    float outFreq = 0.f;
    float peakAmp = 0.f;
    if (binNumber != 0)
    {
        for (int i = 0; i < binNumber; i++)
        {
            int currentBin = lowBin + i;
            float real = out[currentBin][0];
            float imag = out[currentBin][1];

            float mag = std::sqrt(imag * imag + real * real);
            float amp = mag_to_amp(mag, n);

            if (amp > peakAmp)
            {
                outBF.bin = currentBin;
                outBF.freq = binToFreq(currentBin, n, sr);
                peakAmp = amp;
            }

        }
    }
    else
    {
        outFreq = targetFreq;
        outBF.freq = outFreq;
        outBF.bin = freqToBin(outFreq, n, sr);
    }

    return outBF;
}

Sitrano::FreqUnit Sitrano::findPeak(Sitrano::BinFreq inTarget, void* fftwfOut,
    int nfft, float fs, int binRange)
{
    Sitrano::BinFreq target = inTarget;
    Sitrano::FreqUnit outData;

    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(fftwfOut);
    Sitrano::BinFreq outBin;
    outBin.bin = target.bin;

    float log_km1 = logf(hypotf(out[target.bin - 1][0],
        out[target.bin - 1][1]));
    float log_k = logf(hypotf(out[target.bin][0],
        out[target.bin][1]));
    float log_kp1 = logf(hypotf(out[target.bin + 1][0],
        out[target.bin + 1][1]));

    int iStart = -std::abs(binRange);
    int iEnd = std::abs(binRange) + 1;

    if (log_km1 > log_k || log_kp1 > log_k)
    {
        std::vector<float> inputvector;
        for (int i = iStart; i < iEnd; i++)
        {
            inputvector.push_back(logf(hypotf(out[target.bin + i][0],
                out[target.bin + i][1])));
        }
        target.bin += Sitrano::findPeakIndexVector(inputvector) + iStart;

        float log_km1 = logf(hypotf(out[target.bin - 1][0],
            out[target.bin - 1][1]));
        float log_k = logf(hypotf(out[target.bin][0],
            out[target.bin][1]));
        float log_kp1 = logf(hypotf(out[target.bin + 1][0],
            out[target.bin + 1][1]));
    }

    float delta = 0.5f * (log_km1 - log_kp1) / 
        (log_km1 - 2.0f * log_k + log_kp1);

    outBin.fBin = (float)target.bin + delta;
    outBin.freq = binToFreq(outBin.fBin, nfft, fs);

    outData.bin = outBin;

    float interpLogMag = log_k - 0.25f * (log_km1 - log_kp1) * delta;
    outData.mag = expf(interpLogMag);
    outData.amp = mag_to_amp(outData.mag, nfft);

    float phaseAtBin = atan2f(out[target.bin][1], out[target.bin][0]);
    float phase = phaseAtBin - (PI * outData.bin.fBin);

    outData.pha = fmodf(phase + PI, TWO_PI) - PI;

    return outData;
}

void Sitrano::filterVector(std::vector<float>& input, int filterSize, bool zeroPad)
{
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

void Sitrano::normalizeByMaxAbs(std::vector<float>& vec)
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

std::string Sitrano::getRawFilename(std::string& filename)
{
    std::string rawFilename;
    std::filesystem::path path = filename;
    rawFilename = path.filename().replace_extension("").string();

    return rawFilename;
}

void Sitrano::applyWindow(std::vector<float>& frame) 
{
    int size = frame.size();
    for (int i = 0; i < size; ++i) {
        frame[i] *= 0.5f * (1.0f - cosf(TWO_PI * i / (size - 1)));
    }
}

float Sitrano::getPitchFromFilename(const std::string& filename)
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

void Sitrano::saveHarmonicData(
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

int Sitrano::findPeakIndexVector(const std::vector<float>& input)
{
    int outIndex = 0;
    float maxVal = 0.f;

    for (int i = 0; i < input.size(); i++)
    {
        if (input[i] > maxVal) outIndex = i;
        else continue;
    }

    return outIndex;
}

int Sitrano::findPreviousZero(const std::vector<float>& signal, int startSmaple)
{
    int zeroSample = 0;
    for (int i = startSmaple; i > 2 && i < signal.size(); i--)
    {
        bool zeroFound = false;
        zeroFound = (signal[i] >= 0.0 && signal[i - 2] <= 0.0)
            || (signal[i] <= 0.0 && signal[i - 2] >= 0.0);

        if (zeroFound) {
            if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
                zeroSample = i - 2;
            else
                zeroSample = i - 1;
            break;
        }
    }
    return zeroSample;
}

int Sitrano::findNextZero(const std::vector<float>& signal, int startSample) {
    int zeroSample = 0;
    for (int i = startSample; i > 1 && i < signal.size(); i++)
    {
        bool zeroFound = false;
        zeroFound = (signal[i] >= 0.0 && signal[i - 2] < 0.0)
            || (signal[i] <= 0.0 && signal[i - 2] > 0.0);

        if (zeroFound) {
            if (std::abs(signal[i - 2]) <= std::abs(signal[i - 1]))
                zeroSample = i - 2;
            else
                zeroSample = i - 1;
            break;
        }
    }
    return zeroSample;
}

int Sitrano::findNearestZero(const std::vector<float>& signal, int startSample) {
    int prevZero = findPreviousZero(signal, startSample);
    int nextZero = findNextZero(signal, startSample);

    int distancePrev = std::abs(startSample - prevZero);
    int distanceNext = std::abs(startSample - nextZero);

    if (distancePrev < distanceNext)
        return prevZero;
    else
        return nextZero;
}

int Sitrano::findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute) {
    int peakSample = startSample;
    float peakValue = -1.0e30f; // Use a very small number, not 0
    if (useAbsolute) peakValue = 0.f;

    int safeEndSample = std::min((int)signal.size(), endSample);
    int safeStartSample = std::max(0, startSample);

    for (int sample = safeStartSample; sample < safeEndSample; sample++)
    {
        if (useAbsolute) {
            if (std::abs(signal[sample]) > peakValue) {
                peakValue = std::abs(signal[sample]);
                peakSample = sample;
            }
        }
        else {
            if (signal[sample] > peakValue) {
                peakValue = signal[sample];
                peakSample = sample;
            }
        }
    }
    return peakSample;
}

std::vector<int> Sitrano::findZeroCrossings(const std::vector<float>& signal, int startSample) {
    std::vector<int> zeroCrossings;
    for (int i = std::max(1, startSample); i < signal.size(); ++i)
    {
        if ((signal[i - 1] < 0 && signal[i] >= 0) ||
            (signal[i - 1] > 0 && signal[i] <= 0))
        {
            // Find the sample closer to zero
            zeroCrossings.push_back(std::abs(signal[i - 1]) < std::abs(signal[i]) ? (i - 1) : i);
        }
    }
    return zeroCrossings;
}

int Sitrano::findNearestCachedZero(const std::vector<int>& zeroCrossings, int sample) {
    if (zeroCrossings.empty()) return sample; // Safety check

    auto it = std::lower_bound(zeroCrossings.begin(), zeroCrossings.end(), sample);
    if (it == zeroCrossings.end()) return zeroCrossings.back();
    if (it == zeroCrossings.begin()) return *it;

    int after = *it;
    int before = *(it - 1);
    return (std::abs(after - sample) < std::abs(before - sample)) ? after : before;
}

int Sitrano::findAbsPeakIndex(const std::vector<float>& vector) {
    if (vector.empty()) return 0;

    std::vector<float> absTransient(vector.size());
    for (int i = 0; i < absTransient.size(); i++)
    {
        absTransient[i] = std::abs(vector[i]);
    }

    int peakSample = std::distance(absTransient.begin(),
        std::max_element(absTransient.begin(), absTransient.end()));
    return peakSample;
}