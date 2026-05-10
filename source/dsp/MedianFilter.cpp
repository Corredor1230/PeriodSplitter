#include"MedianFilter.h"
#include<algorithm>
#include<stdexcept>

MedianFilter::MedianFilter(Sihat::HPSSSettings set) : nfft(set.nfft), hopSize(set.hopSize), filtSize(set.filtSize), exponent(set.exponent)
{
    initfftw(nfft);
}

void MedianFilter::initfftw(const int nfft)
{
    //time->frequency inits for fft
    inBuffer = (float*)fftwf_malloc(sizeof(float) * nfft);
    outBuffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));

    //frequency->time inits for ifft
    ioutBuffer = (float*)fftwf_malloc(sizeof(float) * nfft);
    invBuffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft / 2 + 1));

    //plans, plans, plans
    ifftPlan = fftwf_plan_dft_c2r_1d(nfft, invBuffer, ioutBuffer, FFTW_ESTIMATE);
    fftPlan = fftwf_plan_dft_r2c_1d(nfft, inBuffer, outBuffer, FFTW_MEASURE);
}

MedianFilter::~MedianFilter()
{
    fftwf_free(inBuffer);
    fftwf_free(outBuffer);
    fftwf_free(ioutBuffer);
    fftwf_free(invBuffer);
    fftwf_destroy_plan(fftPlan);
    fftwf_destroy_plan(ifftPlan);
}

std::vector<float> MedianFilter::processAudio(const std::vector<float>& input)
{

    //first step, getting a regular spectrogram
    Sihat::ComplexSpectrogram spectrogram = getComplexSpectrogram(input);

    //then we filter it (this is the actual median filter)
    MedianFilter::HPSpectrogram filtered = filter(spectrogram);

    //then we apply a mask, basically, we produce two complex spectrograms that have been reduced to vertical or horizontal elements
    MedianFilter::ComplexHPSpec masked = applyMask(filtered, spectrogram);

    //these are the harmonic or percussive audios
    std::vector<float> hAudio = reconstructAudio(masked.harmonic);
    std::vector<float> pAudio = reconstructAudio(masked.percussive);

    //I decided to interleave them for processing simplicity
    std::vector<float> interleavedAudio(hAudio.size() + pAudio.size(), 0.0f);
    for (int i = 0; i < interleavedAudio.size(); i++)
    {
        if (i % 2 == 0) interleavedAudio[i] = hAudio[i / 2];
        else if (i % 2 == 1) interleavedAudio[i] = pAudio[(i - 1) / 2];
    }

    return interleavedAudio;
}

MedianFilter::HPSpectrogram MedianFilter::filter(const Sihat::ComplexSpectrogram& input)
{

    if (filtSize % 2 == 0) 
    {
        std::cerr << "Warning: filtSize MUST be odd. Auto-corrected even number. Evenness is unacceptable here." << std::endl;
        filtSize++;
    }

    //magspec is basically a non-complex spectrogram of magnitudes
    Sihat::Spectrogram magSpec;
    magSpec.numBins = input.numBins;
    magSpec.numFrames = input.numFrames;
    magSpec.data.resize(input.numBins*input.numFrames);

    for (size_t i = 0; i < input.data.size(); ++i) {
        magSpec.data[i] = std::norm(input.data[i]);
    }

    //these are also non-complex
    Sihat::Spectrogram harmonic = magSpec;
    Sihat::Spectrogram percussive = magSpec;

    int radius = filtSize / 2;

    std::vector<float> row(input.numFrames);
    std::vector<float> filteredRow(input.numFrames);
    std::vector<float> hWindow;
    hWindow.reserve(filtSize);

    std::vector<float> pWindow;
    pWindow.resize(filtSize);

    //b stands for bins
    for (int b = 0; b < input.numBins; b++)
    {

        for (int f = 0; f < input.numFrames; f++)
        {
            row[f] = magSpec.at(f, b);
        }

        //f stands for frames, meaning this section handles harmonic stuff
        hWindow.clear();
        for(int i = -radius; i <= radius; i++)
        {
            int idx = std::clamp(i, 0, input.numFrames - 1);
            hWindow.insert(
                std::upper_bound(hWindow.begin(), hWindow.end(), row[idx]),
                row[idx]
            );
        }
        filteredRow[0] = hWindow[radius];

        for(int f = 1; f < input.numFrames; f++)
        {
            int oldIdx = std::clamp(f - 1 - radius, 0, input.numFrames - 1);
            float oldVal = row[oldIdx];

            int newIdx = std::clamp(f + radius, 0, input.numFrames - 1);
            float newVal = row[newIdx];

            if (oldVal != newVal)
            {
                auto it = std::lower_bound(hWindow.begin(), hWindow.end(), oldVal);
                if (it != hWindow.end() && *it == oldVal)
                {
                    hWindow.erase(it);
                }

                hWindow.insert(
                    std::upper_bound(hWindow.begin(), hWindow.end(), newVal),
                    newVal
                );
            }

            filteredRow[f] = hWindow[radius];
        }

        for (int f = 0; f < input.numFrames; f++)
        {
            harmonic.at(f, b) = filteredRow[f];
        }
    }

    //since we're doing frames and then bins, this is vertical percussive stuff
    for (int f = 0; f < input.numFrames; f++)
    {
        for (int b = 0; b < input.numBins; b++)
        {
            for (int i = -radius; i <= radius; i++)
            {
                int idx = std::clamp(b + i, 0, input.numBins - 1);
                pWindow[i + radius] = magSpec.at(f, idx);
            }
            std::nth_element(pWindow.begin(), pWindow.begin() + radius, pWindow.end());
            percussive.at(f, b) = pWindow[radius];
        }
    }

    //HPSpectrogram is just two spectrograms in a single struct. might be useless but I like it because it lets me return stuff.
    return MedianFilter::HPSpectrogram(harmonic, percussive);
}

Sihat::ComplexSpectrogram MedianFilter::getComplexSpectrogram(const std::vector<float>& input)
{
    int numBins = nfft / 2 + 1;
    int numFrames = (input.size() > nfft) ? (input.size() - nfft) / hopSize + 1 : 1;

    Sihat::ComplexSpectrogram spectrogram;
    spectrogram.numFrames = numFrames;
    spectrogram.numBins = numBins;

    spectrogram.data.assign(numFrames * numBins, std::complex<float>(0, 0));

    std::vector<std::complex<float>> frameResult(nfft / 2 + 1);
    for (int frame = 0; frame < numFrames; frame++)
    {
        int frameStart = frame * hopSize;

        for (int samp = 0; samp < nfft; samp++)
        {
            int currSamp = frameStart + samp;
            if (currSamp < input.size())
            inBuffer[samp] = input[currSamp];
            else
            inBuffer[samp] = 0.0;
        }
        Sihat::applyHann(inBuffer, nfft);

        fftwf_execute(fftPlan);

        int frameOffset = frame * numBins;

        for (int i = 0; i <= nfft / 2; i++) {
            spectrogram.data[frameOffset + i] = std::complex<float>(outBuffer[i][0], outBuffer[i][1]);
        }
    }

    return spectrogram;
}

MedianFilter::ComplexHPSpec MedianFilter::applyMask( HPSpectrogram& hp, Sihat::ComplexSpectrogram& complex)
{
    if (hp.harmonic.numFrames != complex.numFrames || hp.percussive.numFrames != complex.numFrames) std::cerr << "All spectrograms should be the same size!!" << std::endl;

    Sihat::ComplexSpectrogram hComplex;
    hComplex.data.resize(complex.data.size());
    hComplex.numBins = complex.numBins;
    hComplex.numFrames = complex.numFrames;
    
    Sihat::ComplexSpectrogram pComplex;
    pComplex.data.resize(complex.data.size());
    pComplex.numBins = complex.numBins;
    pComplex.numFrames = complex.numFrames;

    float hMask = 0.0;
    float pMask = 0.0;
    const float eps = 1e-9f; // Stability constant to prevent NaN

    for (int f = 0; f < complex.numFrames; f++)
    {
        for (int b = 0; b < complex.numBins; b++)
        {
            float hVal = hp.harmonic.at(f, b);
            float pVal = hp.percussive.at(f, b);
            float denominator = std::pow<float>(hVal, exponent) + std::pow<float>(pVal, exponent) + eps;

            float hMask = std::pow<float>(hVal, exponent) / denominator;
            float pMask = std::pow<float>(pVal, exponent) / denominator;

            std::complex<float> original = complex.at(f, b);
            hComplex.at(f, b) = original * hMask;
            pComplex.at(f, b) = original * pMask;
        }
    }

    return ComplexHPSpec(hComplex, pComplex);
}

std::vector<float> MedianFilter::reconstructAudio(const Sihat::ComplexSpectrogram& spec)
{
    int totalSamples = (spec.numFrames - 1) * hopSize + nfft;
    std::vector<float> output(totalSamples, 0.0f);
    std::vector<float> windowSum(totalSamples, 0.0f);

    for (int f = 0; f < spec.numFrames; f++){
        for (int b = 0; b < spec.numBins; b++){
            invBuffer[b][0] = spec.data[f * spec.numBins + b].real();
            invBuffer[b][1] = spec.data[f * spec.numBins + b].imag();
        }

        fftwf_execute(ifftPlan);

        int frameStart = f * hopSize;
        for (int i = 0; i < nfft; i++){
            float windowVal = Sihat::getHannValue(i, nfft);

            float sample = (ioutBuffer[i] / (float)nfft) * windowVal;

            if (frameStart + i < totalSamples) {
                output[frameStart + i] += sample;
                windowSum[frameStart + i] += (windowVal * windowVal);
            }
        }
    }

    for (int i = 0; i < totalSamples; i++){
        if (windowSum[i] > 0.001f){
            output[i] /= windowSum[i];
        } else {
            output[i] = 0.0f;
        }
    }

    return output;
}