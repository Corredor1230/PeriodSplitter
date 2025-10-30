#include "SitranoHeader.h"
#include<iostream>
#include<fftw3.h>
#include<fstream>
#include<filesystem>
#include<iomanip>

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

