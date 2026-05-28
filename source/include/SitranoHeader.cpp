#include "SitranoHeader.h"
#include<iostream>
#include<fftw3.h>
#include<fstream>
#include<filesystem>
#include<iomanip>
#include<complex>

Sihat::BinFreq Sihat::findPeakWithinTolerance(float targetFreq, float tolerance, int n, float sr, void* generic_out, bool isTolInHz = true)
{
    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(generic_out);

    float lowFreq, hiFreq;
    if (!isTolInHz) {
        float midiTolerance = tolerance / 100.0f;
        lowFreq = midiToFreq(freqToMidi(targetFreq) - midiTolerance);        
        hiFreq  = midiToFreq(freqToMidi(targetFreq) + midiTolerance);
    } else {
        lowFreq = targetFreq - tolerance;
        hiFreq  = targetFreq + tolerance;
    }

    // Ensure we don't go out of bounds
    int lowBin = std::max(0, freqToBin(lowFreq, n, sr));
    int hiBin  = std::min(n / 2, freqToBin(hiFreq, n, sr));
    
    // +1 ensures we check the top boundary bin as well
    int binNumber = (hiBin - lowBin) + 1; 

    Sihat::BinFreq outBF{ 0.f, 0.f };
    float peakAmp = -1.0f;
    bool foundPeak = false;

    if (binNumber > 0) {
        for (int i = -1; i < binNumber + 1; i++) {
            int currentBin = lowBin + i;
            if (currentBin < 0) currentBin = 0;
            if (i >= n / 1) break;

            float real = out[currentBin][0];
            float imag = out[currentBin][1];

            float mag = std::sqrt(imag * imag + real * real);
            float amp = mag_to_amp(mag, n);

            // Test begins here
            // std::vector<float> testAmps;
            // for (int b = 0; b < n / 2; b++)
            // {
            //     float testR = out[b][0];
            //     float testI = out[b][1];

            //     float testmag = std::sqrt(testI * testI + testR * testR);
            //     float testamp = mag_to_amp(testmag, n);

            //     testAmps.push_back(testamp);
            // }
            // Test ends here

            if (amp > peakAmp) {
                foundPeak = true;
                outBF.bin = currentBin;
                outBF.freq = binToFreq(currentBin, n, sr);
                peakAmp = amp;
            }
        }
    }
    
    // Fallback if tolerance was too tight to capture a bin
    if (!foundPeak) {
        outBF.freq = targetFreq;
        outBF.bin = freqToBin(targetFreq, n, sr);
    }

    return outBF;
}

Sihat::FreqUnit Sihat::findPeak(Sihat::BinFreq inTarget, void* fftwfOut, int nfft, float fs)
{
    Sihat::FreqUnit outData;
    fftwf_complex* out = reinterpret_cast<fftwf_complex*>(fftwfOut);
    
    int targetBin = inTarget.bin;

    // Safety check against absolute edges
    if (targetBin <= 0 || targetBin >= nfft / 2) {
        outData.bin.bin = targetBin;
        outData.bin.fBin = (float)targetBin;
        outData.bin.freq = binToFreq(targetBin, nfft, fs);
        outData.amp = 0.0f; // Or handle edge case
        return outData;
    }

    float log_km1 = logf(hypotf(out[targetBin - 1][0], out[targetBin - 1][1]));
    float log_k   = logf(hypotf(out[targetBin][0],     out[targetBin][1]));
    float log_kp1 = logf(hypotf(out[targetBin + 1][0], out[targetBin + 1][1]));

    float delta = 0.0f;

    // Only interpolate if we are actually at a local maximum!
    // If not, it means spectral leakage is overpowering this bin. 
    if (log_k > log_km1 && log_k > log_kp1) {
        delta = 0.5f * (log_km1 - log_kp1) / (log_km1 - 2.0f * log_k + log_kp1);
        
        // Clamp delta just in case of float anomalies
        if (delta > 1.0f || delta < -1.0f) delta = 0.0f; 
    }

    outData.bin.bin = targetBin;
    outData.bin.fBin = (float)targetBin + delta;
    outData.bin.freq = binToFreq(outData.bin.fBin, nfft, fs);

    // Calculate magnitude based on whether we interpolated
    if (delta != 0.0f) {
        float interpLogMag = log_k - 0.25f * (log_km1 - log_kp1) * delta;
        outData.mag = expf(interpLogMag);
    } else {
        outData.mag = expf(log_k);
    }
    
    outData.amp = mag_to_amp(outData.mag, nfft);

    float phaseAtBin = atan2f(out[targetBin][1], out[targetBin][0]);
    float phase = phaseAtBin - (PI * outData.bin.fBin);
    outData.pha = fmodf(phase + PI, TWO_PI) - PI;

    return outData;
}

void Sihat::filterVector(std::vector<float>& input, int filterSize, bool zeroPad)
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

void Sihat::normalizeByMaxAbs(std::vector<float>& vec)
{
    if (vec.empty()) return;

    // Find the maximum absolute value
    float maxAbs = 0.0f;
    for (float x : vec) {
        maxAbs = std::max<float>(maxAbs, std::abs(x));
    }

    if (maxAbs == 0.0f) {
        // Avoid division by zero � all elements are zero
        return;
    }

    // Normalize
    for (float& x : vec) {
        x /= maxAbs;
    }
}

std::string Sihat::getRawFilename(const std::string& filename)
{
    std::string rawFilename;
    std::filesystem::path path = filename;
    rawFilename = path.filename().replace_extension("").string();

    return rawFilename;
}

void Sihat::applyWindow(std::vector<float>& frame) 
{
    int size = frame.size();
    for (int i = 0; i < size; ++i) {
        frame[i] *= 0.5f * (1.0f - cosf(TWO_PI * i / (size - 1)));
    }
}

float Sihat::getPitchFromFilename(const std::string& filename)
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

void Sihat::saveHarmonicData(
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

void Sihat::saveHarmonicDataSihat(
    const HarmonicResults& hResults,
    const Sihat::TransientResults& tResults, // <--- NEW ARGUMENT
    float fundamentalFreq,
    const std::string& outputDirectory,
    const std::string& baseFilename,
    const std::string& extension
) {

    const std::vector<uint32_t>& hInd = hResults.finalSamples;
    const std::vector<std::vector<float>>& hAmps = hResults.amps;
    const std::vector<std::vector<float>>& hFreqs = hResults.freqs;
    const float& hRms = hResults.rms;
    const float& tRms = tResults.rms;

    try {
        std::filesystem::create_directories(outputDirectory);

        std::filesystem::path fullPath(outputDirectory);
        fullPath /= baseFilename;
        fullPath.replace_extension(extension);

        std::ofstream f(fullPath, std::ios::binary);
        if (!f) {
            throw std::runtime_error("Failed to open file for writing: " + fullPath.string());
        }

        // Header: Fundamental Frequency
        {
            f.write(reinterpret_cast<const char*>(&fundamentalFreq), sizeof(float));
        }

        // Block 1: Write indices
        {
            uint32_t n = static_cast<uint32_t>(hInd.size());
            f.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
            if (n > 0) {
                f.write(reinterpret_cast<const char*>(hInd.data()), n * sizeof(uint32_t));
            }
        }

        // Block 2: Write amp data
        {
            uint32_t numHarmonics = static_cast<uint32_t>(hAmps.size());
            f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));
            for (const auto& harmonic : hAmps) {
                uint32_t len = static_cast<uint32_t>(harmonic.size());
                f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                if (len > 0) {
                    f.write(reinterpret_cast<const char*>(harmonic.data()), len * sizeof(float));
                }
            }
        }

        // Block 3: Write frequency data
        {
            uint32_t numHarmonics = static_cast<uint32_t>(hFreqs.size());
            f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));

            for (const auto& freqVec : hFreqs) {
                std::vector<ChangePoint> cps;

                if (!freqVec.empty() && !hInd.empty()) {
                    float prev = freqVec.front();
                    cps.push_back({ hInd.front(), prev, prev / fundamentalFreq });

                    for (size_t i = 1; i < freqVec.size() && i < hInd.size(); ++i) {
                        if (std::fabs(freqVec[i] - prev) > 1e-6f) { 
                            cps.push_back({ hInd[i], freqVec[i], freqVec[i] / fundamentalFreq });
                            prev = freqVec[i];
                        }
                    }
                }

                uint32_t len = static_cast<uint32_t>(cps.size());
                f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));

                if (len > 0) {
                    f.write(reinterpret_cast<const char*>(cps.data()), len * sizeof(ChangePoint));
                }
            }
        }

        // Block 4: Write Transient Scalogram Data (NEW)
        {
            uint32_t tStart = static_cast<uint32_t>(tResults.range.initSample);
            uint32_t tEnd   = static_cast<uint32_t>(tResults.range.endSample);
            
            f.write(reinterpret_cast<const char*>(&tStart), sizeof(uint32_t));
            f.write(reinterpret_cast<const char*>(&tEnd), sizeof(uint32_t));

            // B. Write the Scalogram Partials
            uint32_t numPartials = static_cast<uint32_t>(tResults.scalogram.size());
            f.write(reinterpret_cast<const char*>(&numPartials), sizeof(uint32_t));

            for (const auto& partial : tResults.scalogram) {
                // 1. Frequency
                f.write(reinterpret_cast<const char*>(&partial.frequency), sizeof(float));

                // 2. Hop Size 
                // Essential for reconstruction so we know the time delta between envelope points
                uint32_t hop = static_cast<uint32_t>(partial.hopSize);
                f.write(reinterpret_cast<const char*>(&hop), sizeof(uint32_t));

                // 3. Envelope Data Length
                uint32_t dataLen = static_cast<uint32_t>(partial.data.size());
                f.write(reinterpret_cast<const char*>(&dataLen), sizeof(uint32_t));

                // 4. Envelope Data Points
                if (dataLen > 0) {
                    f.write(reinterpret_cast<const char*>(partial.data.data()), dataLen * sizeof(float));
                }
            }
        }

        //Write the RMS values for each section
        {
            f.write(reinterpret_cast<const char*>(&hRms), sizeof(float));
            f.write(reinterpret_cast<const char*>(&tRms), sizeof(float));
        }

        std::cout << "Successfully saved data to: " << fullPath.string() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error saving harmonic data: " << e.what() << std::endl;
    }
}

void Sihat::saveHarmonicDataSihat(
    const HarmonicResults& hResults,
    const Sihat::STransientResults& stResults,
    const Sihat::AnalysisConfig& config, // <--- NEW ARGUMENT
    const Sihat::Settings& aSettings,
    const float fundamentalFreq,
    const uint32_t sr,
    const std::string& outputDirectory,
    const std::string& prefix,
    const std::string& baseFilename,
    const std::string& extension
) {

    const std::vector<uint32_t>& hInd = hResults.finalSamples;
    const std::vector<std::vector<float>>& hAmps = hResults.amps;
    const std::vector<std::vector<float>>& hPhases = hResults.phases;
    const std::vector<std::vector<float>>& hFreqs = hResults.freqs;
    const float& hRms = hResults.rms;
    const float& tRms = stResults.rms;

    try {
        std::filesystem::create_directories(outputDirectory);

        std::filesystem::path fullPath(outputDirectory);
        std::string namewprefix = prefix + baseFilename;
        fullPath /= namewprefix;
        fullPath.replace_extension(extension);

        std::ofstream f(fullPath, std::ios::binary);
        if (!f) {
            throw std::runtime_error("Failed to open file for writing: " + fullPath.string());
        }

        // Header: Sample rate, fundamental frequency and export structure
        {

            //Here's the main structure (did we do transient and harmonic or only one?)
            uint8_t transientAnalysis = aSettings.transientAnalysis ? 1 : 0;
            f.write(reinterpret_cast<const char*>(&transientAnalysis), sizeof(uint8_t));

            uint8_t harmonicAnalysis = aSettings.harmonicAnalysis ? 1 : 0;
            f.write(reinterpret_cast<const char*>(&harmonicAnalysis), sizeof(uint8_t));

            f.write(reinterpret_cast<const char*>(&sr), sizeof(uint32_t));
            f.write(reinterpret_cast<const char*>(&fundamentalFreq), sizeof(float));
            uint32_t filenameLen = static_cast<uint32_t>(baseFilename.size());
            f.write(reinterpret_cast<const char*>(&filenameLen), sizeof(uint32_t));

            if (filenameLen > 0) {f.write(baseFilename.data(), filenameLen);}
            //This is where the export structure would go when I'm finished working it out, maybe a std::vector<bool> but I'm not sure yet.
        }

        if (aSettings.harmonicAnalysis)
        {
            // Block 1: Write indices
            {
                uint32_t n = static_cast<uint32_t>(hInd.size());
                f.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
                if (n > 0) {
                    f.write(reinterpret_cast<const char*>(hInd.data()), n * sizeof(uint32_t));
                }
            }

            // Block 2: Write amp data
            {
                uint32_t numHarmonics = static_cast<uint32_t>(hAmps.size());
                f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));
                for (const auto& harmonic : hAmps) {
                    uint32_t len = static_cast<uint32_t>(harmonic.size());
                    f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    if (len > 0) {
                        f.write(reinterpret_cast<const char*>(harmonic.data()), len * sizeof(float));
                    }
                }
            }

            // Block 2b: Write phase data
            {
                uint32_t numPhases = static_cast<uint32_t>(hPhases.size());
                f.write(reinterpret_cast<const char*>(&numPhases), sizeof(uint32_t));
                for (const auto& phaseList : hPhases) {
                    uint32_t len = static_cast<uint32_t>(phaseList.size());
                    f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    if (len > 0) {
                        f.write(reinterpret_cast<const char*>(phaseList.data()), len * sizeof(float));
                    }
                }
            }

            // Block 3: Write frequency data (ChangePoints)
            {
                uint32_t numHarmonics = static_cast<uint32_t>(hFreqs.size());
                f.write(reinterpret_cast<const char*>(&numHarmonics), sizeof(uint32_t));

                for (const auto& freqVec : hFreqs) {
                    uint32_t len = static_cast<uint32_t>(freqVec.size());
                    f.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
                    if (len > 0){
                        f.write(reinterpret_cast<const char*>(freqVec.data()), len * sizeof(float));
                    }
                }
            }
        }

        if (aSettings.transientAnalysis)
        {
            // Block 4: Write Transient Data
            {
                // 1. Write the Range
                uint32_t tStart = static_cast<uint32_t>(stResults.range.initSample);
                uint32_t tEnd   = static_cast<uint32_t>(stResults.range.endSample);
                
                f.write(reinterpret_cast<const char*>(&tStart), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&tEnd), sizeof(uint32_t));

                // Important parameters (Note: envHopSize is saved right here)
                f.write(reinterpret_cast<const char*>(&stResults.envHopSize), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&stResults.specHopSize), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&stResults.floorHopSize), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&stResults.specWindowSize), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&stResults.specNumBins), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&stResults.specFrameNum), sizeof(uint32_t));
                f.write(reinterpret_cast<const char*>(&config.stSettings.nfft), sizeof(uint32_t));

                // 2. Write scalar metrics
                int32_t riseTime = static_cast<int32_t>(stResults.riseTime);
                f.write(reinterpret_cast<const char*>(&riseTime), sizeof(int32_t));
                f.write(reinterpret_cast<const char*>(&stResults.peakAmp), sizeof(float));

                // Helper lambda for standard float vectors
                auto writeFloatVector = [&f](const std::vector<float>& vec) {
                    uint32_t size = static_cast<uint32_t>(vec.size());
                    f.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
                    if (size > 0) {
                        f.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(float));
                    }
                };

                // 3. Write 1D analysis envelopes
                
                // --- OPTIMIZATION START ---
                // Extract and write only floats for ampEnvelope, saving memory space
                uint32_t envSize = static_cast<uint32_t>(stResults.ampEnvelope.size());
                f.write(reinterpret_cast<const char*>(&envSize), sizeof(uint32_t));
                
                if (envSize > 0) {
                    // Write the initial starting index for perfect reconstruction alongside envHopSize
                    int32_t firstIndex = stResults.ampEnvelope[0].index;
                    f.write(reinterpret_cast<const char*>(&firstIndex), sizeof(int32_t));

                    // Strip out just the float values
                    std::vector<float> envValues(envSize);
                    for (uint32_t i = 0; i < envSize; ++i) {
                        envValues[i] = stResults.ampEnvelope[i].value;
                    }
                    
                    // Write contiguous floats in one I/O operation
                    f.write(reinterpret_cast<const char*>(envValues.data()), envSize * sizeof(float));
                }
                // --- OPTIMIZATION END ---

                writeFloatVector(stResults.centroid);
                writeFloatVector(stResults.flatness);

                // 4. Write the Band Partials
                uint32_t numPartials = static_cast<uint32_t>(config.stSettings.numBands);
                f.write(reinterpret_cast<const char*>(&numPartials), sizeof(uint32_t));
                
                // 5. Write the flattened band envelopes
                writeFloatVector(stResults.bandEnvelopes);

                //This is the metadata for the harmonic section
                uint32_t harmHopSize = static_cast<uint32_t>(stResults.tHarmonics.hopSize);
                f.write(reinterpret_cast<const char*>(&harmHopSize), sizeof(uint32_t));
                uint32_t harmStart = static_cast<uint32_t>(stResults.tHarmonics.startSample);
                f.write(reinterpret_cast<const char*>(&harmStart), sizeof(uint32_t));

                uint32_t numOvertones = static_cast<uint32_t>(stResults.tHarmonics.overtones.size());
                f.write(reinterpret_cast<const char*>(&numOvertones), sizeof(uint32_t));

                //This is the harmonic section proper
                for (const auto& over : stResults.tHarmonics.overtones) {

                    //This is each overtone's target
                    f.write(reinterpret_cast<const char*>(&over.targetOvertone.freq), sizeof(double));
                    f.write(reinterpret_cast<const char*>(&over.targetOvertone.mag), sizeof(double));
                    f.write(reinterpret_cast<const char*>(&over.targetOvertone.amp), sizeof(double));
                    f.write(reinterpret_cast<const char*>(&over.targetOvertone.pha), sizeof(double));

                    uint32_t frameNum = static_cast<uint32_t>(over.envelope.size());
                    f.write(reinterpret_cast<const char*>(&frameNum), sizeof(uint32_t));

                    //This is each overtone's envelope
                    for (const auto& point : over.envelope) {
                        f.write(reinterpret_cast<const char*>(&point.freq), sizeof(float));
                        f.write(reinterpret_cast<const char*>(&point.crestFactor), sizeof(float));
                        f.write(reinterpret_cast<const char*>(&point.amp), sizeof(float));
                    }
                }

                //These are the floor values
                writeFloatVector(stResults.tHarmonics.floor);

                // 6. Write Trajectories
                // uint32_t numTraj = static_cast<uint32_t>(stResults.tHarmonics.overtones.size());
                // f.write(reinterpret_cast<const char*>(&numTraj), sizeof(uint32_t));
                
                // for (const auto& traj : stResults.tHarmonics.overtones) {
                //     f.write(reinterpret_cast<const char*>(&traj.targetOvertone.freq), sizeof(double));
                //     f.write(reinterpret_cast<const char*>(&traj.targetOvertone.mag), sizeof(double));
                //     f.write(reinterpret_cast<const char*>(&traj.targetOvertone.amp), sizeof(double));

                //     f.write(reinterpret_cast<const char*>(&traj.floor), sizeof(float));
                    
                //     uint32_t envLen = static_cast<uint32_t>(traj.envelope.size());
                //     f.write(reinterpret_cast<const char*>(&envLen), sizeof(uint32_t));
                //     if (envLen > 0) {
                //         //f.write(reinterpret_cast<const char*>(traj.envelope.data()), envLen * sizeof(Sihat::TrackedPoint));

                //         for (const auto& point : traj.envelope) {
                //             f.write(reinterpret_cast<const char*>(&point.sampleIndex), sizeof(int32_t));
                //             f.write(reinterpret_cast<const char*>(&point.freq), sizeof(float));
                //             f.write(reinterpret_cast<const char*>(&point.crestFactor), sizeof(float));
                //         }
                //     }
                // }
            }
        }

        // Block 5: Write the RMS values for each section
        {
            if (aSettings.harmonicAnalysis) f.write(reinterpret_cast<const char*>(&hRms), sizeof(float));
            if (aSettings.transientAnalysis) f.write(reinterpret_cast<const char*>(&tRms), sizeof(float));
        }

        std::cout << "Successfully saved data to: " << fullPath.string() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error saving harmonic data: " << e.what() << std::endl;
    }
}

int Sihat::findPeakIndexVector(const std::vector<float>& input)
{
    int outIndex = 0;
    float maxVal = 0.f;

    for (int i = 0; i < input.size(); i++)
    {
        if (input[i] > maxVal)
        {
            maxVal = input[i];
            outIndex = i;
        }
        else continue;
    }

    return outIndex;
}

int Sihat::findPreviousZero(const std::vector<float>& signal, int startSmaple)
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

int Sihat::findNextZero(const std::vector<float>& signal, int startSample) {
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

int Sihat::findNearestZero(const std::vector<float>& signal, int startSample) {
    int prevZero = findPreviousZero(signal, startSample);
    int nextZero = findNextZero(signal, startSample);

    int distancePrev = std::abs(startSample - prevZero);
    int distanceNext = std::abs(startSample - nextZero);

    if (distancePrev < distanceNext)
        return prevZero;
    else
        return nextZero;
}

int Sihat::findPeakSample(const std::vector<float>& signal, int startSample, int endSample, bool useAbsolute) {
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

std::vector<int> Sihat::findZeroCrossings(const std::vector<float>& signal, int startSample) {
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

int Sihat::findNearestCachedZero(const std::vector<int>& zeroCrossings, int sample) {
    if (zeroCrossings.empty()) return sample; // Safety check

    auto it = std::lower_bound(zeroCrossings.begin(), zeroCrossings.end(), sample);
    if (it == zeroCrossings.end()) return zeroCrossings.back();
    if (it == zeroCrossings.begin()) return *it;

    int after = *it;
    int before = *(it - 1);
    return (std::abs(after - sample) < std::abs(before - sample)) ? after : before;
}

int Sihat::findAbsPeakIndex(const std::vector<float>& vector) {
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