#include"FIRFilter.h"

FIRHighPass::FIRHighPass(float sampleRate, float cutoffFreq, float transitionWidth)
{
    designFilter(sampleRate, cutoffFreq, transitionWidth);
    history.assign(taps.size(), 0.0);
}

double FIRHighPass::processSample(float input)
{
    history.insert(history.begin(), input);
    history.pop_back();

    float output = 0.0;
    for (size_t i = 0; i < taps.size(); ++i)
        output += taps[i] * history[i];
    return output;
}

std::vector<float> FIRHighPass::processBuffer(const std::vector<float>& input)
{
    std::vector<float> output(input.size());
    for (size_t i = 0; i < input.size(); ++i)
        output[i] = processSample(input[i]);
    return output;
}

void FIRHighPass::processAudioFile(std::vector<float>& input, const int bufferSize)
{
    int sizeDifference = input.size() % bufferSize;
    for (int i = 0; i < sizeDifference; i++)
    {
        input.push_back(0.0);
    }

    int bufferNumber = input.size() / bufferSize;
    std::vector<float> audioBuffer(bufferSize);
    for (int i = 0; i < bufferNumber; i++)
    {
        for (int samp = 0; samp < bufferSize; samp++)
        {
            int currentSample = i * bufferSize + samp;
            audioBuffer[samp] = input[currentSample];
        }

        audioBuffer = processBuffer(audioBuffer);

        for (int samp = 0; samp < bufferSize; samp++)
        {
            int currentSample = i * bufferSize + samp;
            input[currentSample] = audioBuffer[samp];
        }
    }
}

const std::vector<float>& FIRHighPass::getTaps() const
{
    return taps;
}

void FIRHighPass::designFilter(float sampleRate, float cutoffFreq, float transitionWidth)
{
    double normCutoff = cutoffFreq / sampleRate;
    double normWidth = transitionWidth / sampleRate;

    // Estimate number of taps
    int numTaps = static_cast<int>(std::ceil(4.0 / normWidth));
    if (numTaps % 2 == 0) numTaps++;  // Ensure odd for symmetry

    taps.resize(numTaps);
    int M = numTaps - 1;

    double sum = 0.0;
    for (int n = 0; n <= M; ++n) {
        double centeredN = n - M / 2.0;

        // Ideal lowpass (sinc) then spectral inversion to get highpass
        double sinc = centeredN == 0.0 ?
            2.0 * normCutoff :
            std::sin(2.0 * PI * normCutoff * centeredN) / (PI * centeredN);

        // Hann window
        double window = 0.5 * (1.0 - std::cos(2.0 * PI * n / M));
        taps[n] = sinc * window;
        sum += taps[n];
    }

    for (int n = 0; n <= M; ++n) {
        taps[n] = -taps[n];
    }
    taps[M / 2] += 1.0;
}