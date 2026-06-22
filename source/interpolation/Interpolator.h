#pragma once

#include<vector>
#include<string>

#include"include/InterpHeader.h"
#include"include/InterpConfig.h"
#include"include/ResynthHeader.h"
#include"include/ResynthConfig.h"

class Interpolator
{
public:
    Interpolator(const InterpConfig::Config& config, const InterpConfig::Parameters& p) : conf(config), params(p) {}
    ~Interpolator(){}

    Synth::Sihat interpolateSihat(const Synth::Sihat& a, const Synth::Sihat& b, float f0);

private:
    const InterpConfig::Config& conf;
    const InterpConfig::Parameters& params;
    float dbSilence = -160.0f; // Absolute silence in dB, used for safe padding
    const float ampSilence = 0.0f; // Absolute silence in linear amplitude

    Synth::SihatHeader interpolateHeader(const Synth::SihatHeader& a, const Synth::SihatHeader& b, float f0);
    Synth::SihatHarmonic interpolateHarmonics(const Synth::SihatHarmonic& a, const Synth::SihatHarmonic& b);
    Synth::STransient interpolateTransients(const Synth::STransient& a, const Synth::STransient& b);
    float findDbSilence(const Synth::SihatHarmonic& a, const Synth::SihatHarmonic& b);

};

namespace { // Anonymous namespace for private file-level helpers

    // --- Core Math Helpers ---
    std::vector<float> scaleVector(const std::vector<float>& vec, float scalar) {
        std::vector<float> result = vec;
        for (float& val : result) val *= scalar;
        return result;
    }

    std::vector<Synth::EnvelopePoint> scaleEnvelopeAmp(const std::vector<Synth::EnvelopePoint>& env, float scalar) {
        std::vector<Synth::EnvelopePoint> scaled = env;
        for (auto& pt : scaled) pt.amp *= scalar; // Frequency remains untouched
        return scaled;
    }

    float getMeanFreq(const std::vector<float>& freqs) {
        if (freqs.empty()) return 0.0f;
        float sum = 0.0f;
        for (float f : freqs) sum += f;
        return sum / freqs.size();
    }

    // --- Matching Struct ---
    struct OvertoneMatch {
        int indexA = -1;
        int indexB = -1;
        bool isMatched = false;
    };

    // --- Harmonic Matcher ---
    std::vector<OvertoneMatch> matchHarmonics(
        const std::vector<std::vector<float>>& freqsA,
        const std::vector<std::vector<float>>& freqsB,
        double tolerance = 0.1) 
    {
        std::vector<OvertoneMatch> matches;
        std::vector<bool> b_used(freqsB.size(), false);

        for (int i = 0; i < freqsA.size(); ++i) {
            float meanA = getMeanFreq(freqsA[i]);
            int best_b_idx = -1;
            double min_dist = tolerance;

            for (int j = 0; j < freqsB.size(); ++j) {
                if (b_used[j]) continue;
                double dist = std::abs(meanA - getMeanFreq(freqsB[j]));
                if (dist < min_dist) {
                    min_dist = dist;
                    best_b_idx = j;
                }
            }

            if (best_b_idx != -1) {
                matches.push_back({i, best_b_idx, true});
                b_used[best_b_idx] = true;
            } else {
                matches.push_back({i, -1, false});
            }
        }

        for (int j = 0; j < freqsB.size(); ++j) {
            if (!b_used[j]) matches.push_back({-1, j, false});
        }
        return matches;
    }

    // --- Transient Matcher ---
    std::vector<OvertoneMatch> matchOvertones(
        const std::vector<Synth::TransientOvertone>& overtonesA,
        const std::vector<Synth::TransientOvertone>& overtonesB,
        double tolerance = 0.1) 
    {
        std::vector<OvertoneMatch> matches;
        std::vector<bool> b_used(overtonesB.size(), false);

        for (int i = 0; i < overtonesA.size(); ++i) {
            int best_b_idx = -1;
            double min_dist = tolerance;

            for (int j = 0; j < overtonesB.size(); ++j) {
                if (b_used[j]) continue;
                double dist = std::abs(overtonesA[i].target.freq - overtonesB[j].target.freq);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_b_idx = j;
                }
            }

            if (best_b_idx != -1) {
                matches.push_back({i, best_b_idx, true});
                b_used[best_b_idx] = true;
            } else {
                matches.push_back({i, -1, false});
            }
        }

        for (int j = 0; j < overtonesB.size(); ++j) {
            if (!b_used[j]) matches.push_back({-1, j, false});
        }
        return matches;
    }


    std::vector<Synth::ModalComponent> matchComponents(
        const std::vector<Synth::ModalComponent>& modesA,
        const std::vector<Synth::ModalComponent>& modesB,
        double tolerance = 10.0) 
    {
        std::vector<Synth::ModalComponent> matched;
        std::vector<bool> b_used(modesB.size(), false);

        for (int i = 0; i < modesA.size(); ++i) {
            int best_b_idx = -1;
            double min_dist = tolerance;

            for (int j = 0; j < modesB.size(); ++j) {
                if (b_used[j]) continue;
                double dist = std::abs(modesA[i].freq - modesB[j].freq);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_b_idx = j;
                }
            }

            if (best_b_idx != -1) {
                matched.push_back(modesA[i]); // Or some combination logic
                b_used[best_b_idx] = true;
            } else {
                matched.push_back(modesA[i]);
            }
        }

        for (int j = 0; j < modesB.size(); ++j) {
            if (!b_used[j]) matched.push_back(modesB[j]);
        }
        return matched;
    }

    template <typename T>
    std::vector<T> padWithValue(const std::vector<T>& vec, size_t targetSize, T padVal) {
        if (vec.size() >= targetSize) return vec;
        std::vector<T> out = vec;
        out.resize(targetSize, padVal);
        return out;
    }

    // Helper: Pads a vector to a target size by repeating its last element (e.g., for frequencies)
    template <typename T>
    std::vector<T> padWithLast(const std::vector<T>& vec, size_t targetSize) {
        if (vec.size() >= targetSize || vec.empty()) return vec;
        std::vector<T> out = vec;
        out.resize(targetSize, vec.back());
        return out;
    }

    struct ModeMatch {
        int indexA = -1;
        int indexB = -1;
        bool isMatched = false;
    };

    std::vector<ModeMatch> matchModalComponents(
        const std::vector<Synth::ModalComponent>& modesA,
        const std::vector<Synth::ModalComponent>& modesB,
        float freqToleranceRatio = 0.05f) // 5% difference threshold (~1 semitone)
    {
        std::vector<ModeMatch> matches;
        std::vector<bool> b_used(modesB.size(), false);

        for (int i = 0; i < modesA.size(); ++i) {
            int best_b_idx = -1;
            // The allowed deviation scales with the frequency
            float min_dist = modesA[i].freq * freqToleranceRatio; 

            for (int j = 0; j < modesB.size(); ++j) {
                if (b_used[j]) continue;
                
                float dist = std::abs(modesA[i].freq - modesB[j].freq);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_b_idx = j;
                }
            }

            if (best_b_idx != -1) {
                matches.push_back({i, best_b_idx, true});
                b_used[best_b_idx] = true;
            } else {
                matches.push_back({i, -1, false});
            }
        }

        for (int j = 0; j < modesB.size(); ++j) {
            if (!b_used[j]) {
                matches.push_back({-1, j, false});
            }
        }
        return matches;
    }

}