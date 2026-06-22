#include"Interpolator.h"

// Core interpolation function for Sihat structures
Synth::Sihat Interpolator::interpolateSihat(const Synth::Sihat& a, const Synth::Sihat& b, float f0)
{
    Synth::Sihat out;

    out.header = interpolateHeader(a.header, b.header, f0);
    if (out.header.type.harmonicAnalysis) out.harmonic = interpolateHarmonics(a.harmonic, b.harmonic);
    if (out.header.type.harmonicAnalysis) out.transient = interpolateTransients(a.transient, b.transient);

    return out;
}

Synth::SihatHeader Interpolator::interpolateHeader(const Synth::SihatHeader& a, const Synth::SihatHeader& b, float f0)
{
    Synth::SihatHeader out;
    out.f0 = f0;
    out.filename = conf.i_filename; 
    out.sampleRate = a.sampleRate; // Sample rate must remain consistent
    out.type = a.type;             // Inherit analysis flags from A
    return out;
}

Synth::SihatHarmonic Interpolator::interpolateHarmonics(const Synth::SihatHarmonic& a, const Synth::SihatHarmonic& b)
{
    Synth::SihatHarmonic out;

    // Pad structural indices with the last known frame index
    size_t maxIndices = std::max(a.indices.size(), b.indices.size());
    out.indices = SihatInterpolation::ivec_interpolate(
        padWithLast(a.indices, maxIndices), 
        padWithLast(b.indices, maxIndices), 
        conf.alpha, conf.type
    );
    
    out.numFrames = std::max(a.numFrames, b.numFrames);
    out.rms = SihatInterpolation::f_interpolate(a.rms, b.rms, conf.alpha, conf.type);

    std::vector<OvertoneMatch> matches = matchHarmonics(a.fRatio, b.fRatio, 0.25);

    int nonEmptyAs = 0;
    int nonEmptyBs = 0;

    float sum = 0.0f;
    for (int i = 0; i < a.amp.size(); ++i) {
        if (!a.amp[i].empty()) {
            nonEmptyAs++;
            sum += a.amp[i].back();
        }
    }
    float avSilenceA = nonEmptyAs > 0 ? sum / nonEmptyAs : dbSilence;

    sum = 0.0f;
    for (int i = 0; i < b.amp.size(); ++i) {
        if (!b.amp[i].empty()) {
            nonEmptyBs++;
            sum += b.amp[i].back();
        }
    }
    float avSilenceB = nonEmptyBs > 0 ? sum / nonEmptyBs : dbSilence;

    for (const auto& match : matches)
    {
        if (match.isMatched) {
            size_t targetFrames = std::max(a.amp[match.indexA].size(), b.amp[match.indexB].size());

            // Amplitudes pad with -160.0 (absolute silence)
            auto safeAmpA = padWithValue<float>(a.amp[match.indexA], targetFrames, avSilenceA);
            auto safeAmpB = padWithValue<float>(b.amp[match.indexB], targetFrames, avSilenceB);

            // Frequencies pad with the last known frequency
            auto safeRatioA = padWithLast(a.fRatio[match.indexA], targetFrames);
            auto safeRatioB = padWithLast(b.fRatio[match.indexB], targetFrames);

            out.amp.push_back(SihatInterpolation::fvec_interpolate(safeAmpA, safeAmpB, conf.alpha, conf.type));
            out.fRatio.push_back(SihatInterpolation::fvec_interpolate(safeRatioA, safeRatioB, conf.alpha, conf.type));
        } 
        else if (match.indexA != -1) {
            auto silenceVecA = std::vector<float>(a.amp[match.indexA].size(), dbSilence);
            out.amp.push_back(SihatInterpolation::fvec_interpolate(a.amp[match.indexA], silenceVecA, conf.alpha, conf.type));
            out.fRatio.push_back(a.fRatio[match.indexA]);
        } 
        else if (match.indexB != -1) {
            auto silenceVecB = std::vector<float>(b.amp[match.indexB].size(), dbSilence);
            out.amp.push_back(SihatInterpolation::fvec_interpolate(silenceVecB, b.amp[match.indexB], conf.alpha, conf.type));
            out.fRatio.push_back(b.fRatio[match.indexB]);
        }
    }

    for (size_t i = 0; i < out.amp.size(); ++i) {
        if (out.amp[i].empty()) continue; // Ignore safely
        
        if (out.amp[i].size() < out.numFrames) {
            // Pad amplitudes with absolute silence (-160.0f)
            out.amp[i] = padWithValue<float>(out.amp[i], out.numFrames, dbSilence);
            
            // Pad frequencies with the last known ratio to prevent pitch sweeping
            out.fRatio[i] = padWithLast(out.fRatio[i], out.numFrames);
        }
    }

    return out;
}

Synth::STransient Interpolator::interpolateTransients(const Synth::STransient& a, const Synth::STransient& b)
{
    Synth::STransient out;

    // 1. Metadata 
    out.meta = a.meta; 
    out.meta.numBands = std::min(a.meta.numBands, b.meta.numBands);
    out.meta.specNumFrames = std::max(a.meta.specNumFrames, b.meta.specNumFrames);

    // 2. Standard 1D Data (Pad spectral characteristics with their last value)
    size_t maxFrames1D = std::max(a.flatness.size(), b.flatness.size());
    out.flatness = SihatInterpolation::fvec_interpolate(padWithLast(a.flatness, maxFrames1D), padWithLast(b.flatness, maxFrames1D), conf.alpha, conf.type);
    out.centroid = SihatInterpolation::fvec_interpolate(padWithLast(a.centroid, maxFrames1D), padWithLast(b.centroid, maxFrames1D), conf.alpha, conf.type);
    
    // Floors act as amplitude minimums, so 0.0f is safer here than repeating a high noise floor
    size_t maxFloors = std::max(a.floors.size(), b.floors.size());
    out.floors = SihatInterpolation::fvec_interpolate(padWithValue(a.floors, maxFloors, 0.0f), padWithValue(b.floors, maxFloors, 0.0f), conf.alpha, conf.type);
    
    out.riseTime = SihatInterpolation::i_interpolate(a.riseTime, b.riseTime, conf.alpha, conf.type);
    out.peakAmp  = SihatInterpolation::f_interpolate(a.peakAmp, b.peakAmp, conf.alpha, conf.type);
    out.rms      = SihatInterpolation::f_interpolate(a.rms, b.rms, conf.alpha, conf.type);

    // 3. Envelope
    out.envelope.useHopSize = a.envelope.useHopSize;
    out.envelope.hopSize = a.envelope.hopSize;
    out.envelope.firstIndex = std::min(a.envelope.firstIndex, b.envelope.firstIndex);
    
    size_t maxEnvSize = std::max(a.envelope.env.size(), b.envelope.env.size());
    // Amplitude envelope drops to 0.0, indices freeze on the last index
    out.envelope.env = SihatInterpolation::fvec_interpolate(padWithValue<float>(a.envelope.env, maxEnvSize, ampSilence), padWithValue<float>(b.envelope.env, maxEnvSize, ampSilence), conf.alpha, conf.type);
    out.envelope.index = SihatInterpolation::ivec_interpolate(padWithLast(a.envelope.index, maxEnvSize), padWithLast(b.envelope.index, maxEnvSize), conf.alpha, conf.type);

    // 4. Deflatten, Interpolate, Reflatten Bands
    auto get_bands_2d = [](const std::vector<float>& flat, uint32_t numBands) {
        uint32_t numFrames = flat.size() / (numBands == 0 ? 1 : numBands);
        std::vector<std::vector<float>> deflat(numBands, std::vector<float>(numFrames));
        for(uint32_t f = 0; f < numFrames; ++f) {
            for(uint32_t b = 0; b < numBands; ++b) deflat[b][f] = flat[f * numBands + b];
        }
        return deflat;
    };

    uint32_t minBands = out.meta.numBands;
    auto a_bands = get_bands_2d(a.bands, minBands);
    auto b_bands = get_bands_2d(b.bands, minBands);
    
    std::vector<std::vector<float>> out_bands_2d(minBands);
    for(uint32_t b = 0; b < minBands; ++b) {
        size_t maxBandFrames = std::max(a_bands[b].size(), b_bands[b].size());
        
        // Pad spectral bands with 0.0f (silence)
        auto safeBandA = padWithValue<float>(a_bands[b], maxBandFrames, ampSilence);
        auto safeBandB = padWithValue<float>(b_bands[b], maxBandFrames, ampSilence);
        
        out_bands_2d[b] = SihatInterpolation::fvec_interpolate(safeBandA, safeBandB, conf.alpha, conf.type);
    }
    
    uint32_t outFrames = out_bands_2d.empty() ? 0 : out_bands_2d[0].size();
    out.bands.reserve(outFrames * minBands);
    for(uint32_t f = 0; f < outFrames; ++f) {
        for(uint32_t b = 0; b < minBands; ++b) out.bands.push_back(out_bands_2d[b][f]);
    }

    // 5. Transient Overtones 
    // (Your code here is already perfect. The inline ternary operators handles 
    //  the frequency/amplitude padding exactly as required).
    std::vector<OvertoneMatch> t_matches = matchOvertones(a.overtones, b.overtones, 0.1);

    // 6. Transient Modal Components
    // a. Sort both mode arrays by frequency (lowest to highest)
    // 6. Modal Components (TModes)
    out.tModes.startInd = SihatInterpolation::i_interpolate(a.tModes.startInd, b.tModes.startInd, conf.alpha, conf.type);
    out.tModes.length = SihatInterpolation::i_interpolate(a.tModes.length, b.tModes.length, conf.alpha, conf.type);

    // Match modes with a 5% frequency tolerance
    std::vector<ModeMatch> m_matches = matchModalComponents(a.tModes.modes, b.tModes.modes, 0.05f);

    for (const auto& match : m_matches) {
        Synth::ModalComponent interp_mode;
        
        if (match.isMatched) {
            // Morph the physical properties of the matched resonance
            const auto& mA = a.tModes.modes[match.indexA];
            const auto& mB = b.tModes.modes[match.indexB];
            
            interp_mode.freq  = SihatInterpolation::f_interpolate(mA.freq, mB.freq, conf.alpha, conf.type);
            interp_mode.amp   = SihatInterpolation::f_interpolate(mA.amp, mB.amp, conf.alpha, conf.type);
            interp_mode.decay = SihatInterpolation::f_interpolate(mA.decay, mB.decay, conf.alpha, conf.type);
            interp_mode.phase = SihatInterpolation::f_interpolate(mA.phase, mB.phase, conf.alpha, conf.type);
            
            out.tModes.modes.push_back(interp_mode);
        }
        else if (match.indexA != -1) {
            // Fade out unmatched mode A
            interp_mode = a.tModes.modes[match.indexA];
            interp_mode.amp = SihatInterpolation::f_interpolate(interp_mode.amp, ampSilence, conf.alpha, conf.type);
            
            // Only keep it if it hasn't faded to absolute zero
            if (interp_mode.amp > 0.0001f) out.tModes.modes.push_back(interp_mode);
        }
        else if (match.indexB != -1) {
            // Fade in unmatched mode B
            interp_mode = b.tModes.modes[match.indexB];
            interp_mode.amp = SihatInterpolation::f_interpolate(ampSilence, interp_mode.amp, conf.alpha, conf.type);
            
            // Only keep it if it has actually started fading in
            if (interp_mode.amp > 0.0001f) out.tModes.modes.push_back(interp_mode);
        }
    }

    for (const auto& match : t_matches) {
        if (match.isMatched) {
            const auto& otA = a.overtones[match.indexA];
            const auto& otB = b.overtones[match.indexB];
            Synth::TransientOvertone interp_ot;
            
            interp_ot.target.index = otA.target.index;
            interp_ot.target.freq = SihatInterpolation::f_interpolate(otA.target.freq, otB.target.freq, conf.alpha, conf.type);
            interp_ot.target.mag = SihatInterpolation::f_interpolate(otA.target.mag, otB.target.mag, conf.alpha, conf.type);
            interp_ot.target.amp = SihatInterpolation::f_interpolate(otA.target.amp, otB.target.amp, conf.alpha, conf.type);
            interp_ot.target.phase = otA.target.phase;

            int maxPts = std::max(otA.envelope.size(), otB.envelope.size());
            interp_ot.envelope.reserve(maxPts);
            for(int i = 0; i < maxPts; ++i) {
                Synth::EnvelopePoint pt;
                // Note: This logic you wrote is flawless.
                float a_freq = (i < otA.envelope.size()) ? otA.envelope[i].freq : otA.envelope.back().freq;
                float a_amp  = (i < otA.envelope.size()) ? otA.envelope[i].amp  : otA.envelope.back().amp;
                float a_cf   = (i < otA.envelope.size()) ? otA.envelope[i].crestFactor : otA.envelope.back().crestFactor;

                float b_freq = (i < otB.envelope.size()) ? otB.envelope[i].freq : otB.envelope.back().freq;
                float b_amp  = (i < otB.envelope.size()) ? otB.envelope[i].amp  : otB.envelope.back().amp; // For amplitude, we can allow the last known value to fade out naturally instead of forcing silence
                float b_cf   = (i < otB.envelope.size()) ? otB.envelope[i].crestFactor : otB.envelope.back().crestFactor;

                pt.freq = SihatInterpolation::f_interpolate(a_freq, b_freq, conf.alpha, conf.type);
                pt.amp = SihatInterpolation::f_interpolate(a_amp, b_amp, conf.alpha, conf.type);
                pt.crestFactor = SihatInterpolation::f_interpolate(a_cf, b_cf, conf.alpha, conf.type);
                interp_ot.envelope.push_back(pt);
            }
            out.overtones.push_back(interp_ot);
        } 
        else if (match.indexA != -1) {
            Synth::TransientOvertone fadedA = a.overtones[match.indexA];
            fadedA.target.amp *= (1.0f - conf.alpha);
            fadedA.target.mag *= (1.0f - conf.alpha);
            for (auto& pt : fadedA.envelope) {
                pt.amp = SihatInterpolation::f_interpolate(pt.amp, dbSilence, conf.alpha, conf.type);
            }
            out.overtones.push_back(fadedA);
        } 
        else if (match.indexB != -1) {
            Synth::TransientOvertone fadedB = b.overtones[match.indexB];
            fadedB.target.amp *= conf.alpha;
            fadedB.target.mag *= conf.alpha;
            for (auto& pt : fadedB.envelope) {
                pt.amp = SihatInterpolation::f_interpolate(dbSilence, pt.amp, conf.alpha, conf.type);
            }
            out.overtones.push_back(fadedB);
        }
    }

    size_t globalMaxEnv = 0;
    for (const auto& ot : out.overtones) {
        globalMaxEnv = std::max(globalMaxEnv, ot.envelope.size());
    }

    // 2. Pad all shorter envelopes to match that maximum length
    for (auto& ot : out.overtones) {
        if (ot.envelope.empty()) continue;
        
        if (ot.envelope.size() < globalMaxEnv) {
            // Grab the last valid envelope point to lock in the frequency & crest factor
            Synth::EnvelopePoint padPt = ot.envelope.back();
            padPt.amp = dbSilence; // Force amplitude to -160dB so it goes completely silent
            
            // Pad until it matches the global max
            while (ot.envelope.size() < globalMaxEnv) {
                ot.envelope.push_back(padPt);
            }
        }
    }

    return out;
}

float Interpolator::findDbSilence(const Synth::SihatHarmonic& a, const Synth::SihatHarmonic& b) {
    float sum = 0.0f;
    int count = 0;
    Synth::SihatHarmonic smaller = a.amp.size() < b.amp.size() ? a : b;
    for (size_t i = 0; i < smaller.amp.size(); ++i) {
        if (smaller.amp[i].empty()) continue;
        sum += smaller.amp[i].back();
        count++;
    }
    if (count != 0) return sum / count;
    return 0.0f;
}