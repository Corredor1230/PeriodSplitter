#pragma once

#include<cstdint>
#include<vector>
#include<string>

namespace Synth
{
    constexpr double PI = 3.14159265358979323846;
    constexpr double M_E = 2.718281828459045;
    constexpr double TWO_PI = 2.0 * PI;

    constexpr char SIHAT_MAGIC[4] = {'S', 'I', 'H', 'T'};
    constexpr uint16_t SIHAT_VERSION_MAJOR = 1;
    constexpr uint16_t SIHAT_VERSION_MINOR = 0;

    struct AnalysisType{
        bool transientAnalysis = true;
        bool harmonicAnalysis = true;
    };

    struct ModalComponent{
        float freq = 0.0f;
        float amp = 0.0f;
        float phase = 0.0f;
        float decay = 0.0f;
    };

    struct TModes{
        std::vector<ModalComponent> modes;
        uint32_t startInd = 0;
        uint32_t length = 1024;
    };

    struct Envelope{
        bool useHopSize = true;
        uint32_t hopSize = 128;
        uint32_t firstIndex = 0;
        std::vector<float> env;
        std::vector<int> index;
    };

    struct SpectralBin{
        uint32_t index = 0;
        double freq = 0.0;
        double mag = 0.0;
        double amp = 0.0;
        float phase = 0.0;
    };

    struct EnvelopePoint{
        float freq = 0.0;
        float amp = 0.0;
        float crestFactor = 0.0;
    };

    struct TransientMetadata{
        uint32_t tStart = 0;
        uint32_t tEnd = 0;
        uint32_t envHopSize = 64;
        uint32_t specHopSize = 64;
        uint32_t floorHopSize = 64;
        uint32_t specWindowSize = 1024;
        uint32_t specNumBins = 16;
        uint32_t specNumFrames = 64;
        uint32_t specNfft = 2048;
        uint32_t numBands = 16;
        uint32_t harmHopSize = 16;
        uint32_t harmStartSample = 1000;
    };

    struct HarmonicMetadata{
        uint32_t hopSize = 32;
        uint32_t startSample = 1000;
    };

    struct TransientOvertone{
        //Goal frequency
        SpectralBin target;
        std::vector<EnvelopePoint> envelope;
    };

    struct SihatHeader{
        uint32_t sampleRate = 96000;
        float f0 = 60.0;
        std::string filename = "data";
        AnalysisType type;
    };

    struct SihatHarmonic{
        uint32_t numFrames = 0;
        std::vector<uint32_t> indices;
        std::vector<std::vector<float>> amp;
        std::vector<std::vector<float>> pha;
        std::vector<std::vector<float>> fRatio;
        float rms;
    };

    struct STransient{
        TransientMetadata meta;
        Envelope envelope;
        std::vector<TransientOvertone> overtones;
        std::vector<float> floors;
        std::vector<float> flatness;
        std::vector<float> centroid;
        std::vector<float> bands;
        TModes tModes;
        int32_t riseTime = 1000;
        float peakAmp = 0.5;
        float rms;
    };

    struct Sihat{
        SihatHeader header;
        SihatHarmonic harmonic;
        STransient transient;
    };

    
}