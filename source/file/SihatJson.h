#pragma once

#include"include/json.hpp"
#include"include/SitranoHeader.h"
#include<fstream>

using json = nlohmann::json;

class SihatJson {
public:
    // ---- PitchSettings ----
    void to_json(json& j, const Sitrano::PitchSettings& p) {
        j = json{
            {"p_threshold", p.modeThreshold},
            {"p_inCents", p.toleranceInCents},
            {"p_minFreq", p.minFreq},
            {"p_maxFreq", p.maxFreq}
        };
    }

    // ---- TransientSettings ----
    void to_json(json& j, const Sitrano::TransientSettings& t) {
        j = json{
            {"t_StartSample", t.tStartSample},
            {"t_UseMs", t.useMs},
            {"t_rmsSampleSize", t.rmsSampleSize},
            {"t_rmsHopLength", t.rmsSampleHopLength},
            {"t_rmsSizeMs", t.transientRmsSizeMs},
            {"t_rmsHopRatio", t.transientRmsHopRatio},
            {"t_factor", t.transientFactor},
            {"t_threshold", t.transientThreshold},
            {"t_preAttack", t.preAttack}
        };
    }

    // ---- CorrelationSettings ----
    void to_json(json& j, const Sitrano::CorrelationSettings& c) {
        j = json{
            {"c_offset", c.periodStartOffsetMs},
            {"c_threshold", c.correlationThreshold}
        };
    }

    // ---- OvertoneSettings ----
    void to_json(json& j, const Sitrano::OvertoneSettings& o) {
        j = json{
            {"o_useTolerance", o.useTolerance},
            {"o_tolerance", o.toleranceValue},
            {"o_postTransientStart", o.postTransientStart},
            {"o_FirstSample", o.overtoneFirstSample},
            {"o_useCustomSignal", o.useCustomSignal},
            {"o_sumAmplitudes", o.sumAmplitudes},
            {"o_mfft", o.fftSize},
            {"o_ignoreThreshold", o.overtoneThreshold},
            {"o_absThreshold", o.setAbsoluteThreshold}
        };
    }

    // ---- HarmonicSettings ----
    void to_json(json& j, const Sitrano::HarmonicSettings& h) {
        j = json{
            {"h_enabled", h.applyHanning},
            {"h_windowStyle", static_cast<int>(h.style)},  // convert enum to int
            {"h_Tolerance", h.toleranceValue}
        };
    }

    // ---- General Settings ----
    void to_json(json& j, const Sitrano::Settings& s) {
        j = json{
            {"pitchAnalysis", s.pitchAnalysis},
            {"transientSeparation", s.transientSeparation},
            {"periodAnalysis", s.periodAnalysis},
            {"overtoneAnalysis", s.overtoneAnalysis},
            {"harmonicAnalysis", s.harmonicAnalysis},
            {"noiseAnalysis", s.noiseAnalysis}
        };
    }

    // ---- AnalysisConfig ----
    void to_json(json& j, const Sitrano::AnalysisConfig& c) {
        j = json{
            {"maxHarmonics", c.numHarmonics},
            {"N", c.nfft},
            {"hopSize", c.hopSize},
            {"startSample", c.startSample},
            {"tolerance", c.tolerance},
            {"pitchSettings", c.pConfig},
            {"transientSettings", c.tSettings},
            {"correlationSettings", c.cSettings},
            {"overtoneSettings", c.oConfig},
            {"harmonicSettings", c.hConfig},
            {"verbose", c.verbose}
        };
    }

};