#pragma once

#include"include/json.hpp"
#include"include/SitranoHeader.h"
#include"file/SihatFile.h"
#include<fstream>

using json = nlohmann::json;

namespace SihatJson {
    // ---- PitchSettings ----
    inline void to_json(json& j, const Sitrano::PitchSettings& p) {
        j = json{
            {"p_threshold", p.modeThreshold},
            {"p_inCents", p.toleranceInCents},
            {"p_minFreq", p.minFreq},
            {"p_maxFreq", p.maxFreq}
        };
    }

    // ---- TransientSettings ----
    inline void to_json(json& j, const Sitrano::TransientSettings& t) {
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

    inline void to_json(json& j, const SihatFile::OutInfo& i) {
        j = json{
            {"outDir", i.outDir},
            {"prefix", i.prefix},
            {"extension", i.extension}
        };
    }

    // ---- CorrelationSettings ----
    inline void to_json(json& j, const Sitrano::CorrelationSettings& c) {
        j = json{
            {"c_offset", c.periodStartOffsetMs},
            {"c_threshold", c.correlationThreshold}
        };
    }

    // ---- OvertoneSettings ----
    inline void to_json(json& j, const Sitrano::OvertoneSettings& o) {
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
    inline void to_json(json& j, const Sitrano::HarmonicSettings& h) {
        j = json{
            {"h_enabled", h.applyHanning},
            {"h_windowStyle", static_cast<int>(h.style)},  // convert enum to int
            {"h_Tolerance", h.toleranceValue}
        };
    }

    // --- NoiseSettings ---
    inline void to_json(json& j, const Sitrano::NoiseSettings& n) {
        j = json{
            {"n_nfft", n.nfft},
            {"n_hopSize", n.hopSize},
            {"n_numBins", n.numBins},
            {"n_minFreq", n.minFreq},
            {"n_maxFreq", n.maxFreq},
            {"n_useOctaveDiv", n.useOctaveDiv},
            {"n_octaveDiv", n.octaveDiv},
            {"n_startSample", n.startSample},
            {"n_useList", n.useList}
        };
    }

    // ---- General Settings ----
    inline void to_json(json& j, const Sitrano::Settings& s) {
        j = json{
            {"fullHarmonicAnalysis", s.fullHarmonicAnalysis},
            {"pitchAnalysis", s.pitchAnalysis},
            {"transientSeparation", s.transientSeparation},
            {"periodAnalysis", s.periodAnalysis},
            {"overtoneAnalysis", s.overtoneAnalysis},
            {"harmonicAnalysis", s.harmonicAnalysis},
            {"noiseAnalysis", s.noiseAnalysis}
        };
    }

    // ---- AnalysisConfig ----
    inline void to_json(json& j, const Sitrano::AnalysisConfig& c) {

        json pSettings;
        json tSettings;
        json cSettings;
        json oSettings;
        json hSettings;
        json nSettings;

        to_json(pSettings, c.pSettings);
        to_json(tSettings, c.tSettings);
        to_json(cSettings, c.cSettings);
        to_json(oSettings, c.oSettings);
        to_json(hSettings, c.hSettings);
        to_json(nSettings, c.nSettings);


        j = json{
            {"maxHarmonics", c.numHarmonics},
            {"N", c.nfft},
            {"hopSize", c.hopSize},
            {"startSample", c.startSample},
            {"tolerance", c.tolerance},
            {"pitchSettings", pSettings},
            {"transientSettings", tSettings},
            {"correlationSettings", cSettings},
            {"overtoneSettings", oSettings},
            {"harmonicSettings", hSettings},
            {"noiseSettings", nSettings},
            {"verbose", c.verbose},
            {"bulkProcess", c.bulkProcess}
        };
    }

    // --- PitchSettings ---
    inline void from_json(const json& j, Sitrano::PitchSettings& p) {
        j.at("p_threshold").get_to(p.modeThreshold);
        j.at("p_inCents").get_to(p.toleranceInCents);
        j.at("p_minFreq").get_to(p.minFreq);
        j.at("p_maxFreq").get_to(p.maxFreq);
    }

    // ---- TransientSettings ----
    inline void from_json(const json& j, Sitrano::TransientSettings& t) {
        j.at("t_StartSample").get_to(t.tStartSample);
        j.at("t_UseMs").get_to(t.useMs);
        j.at("t_rmsSampleSize").get_to(t.rmsSampleSize);
        j.at("t_rmsHopLength").get_to(t.rmsSampleHopLength);
        j.at("t_rmsSizeMs").get_to(t.transientRmsSizeMs);
        j.at("t_rmsHopRatio").get_to(t.transientRmsHopRatio);
        j.at("t_factor").get_to(t.transientFactor);
        j.at("t_threshold").get_to(t.transientThreshold);
        j.at("t_preAttack").get_to(t.preAttack);
    }

    // --- NoiseSettings ---
    inline void from_json(const json& j, Sitrano::NoiseSettings& n) {
        j.at("n_nfft").get_to(n.nfft);
        j.at("n_hopSize").get_to(n.hopSize);
        j.at("n_numBins").get_to(n.numBins);
        j.at("n_minFreq").get_to(n.minFreq);
        j.at("n_maxfreq").get_to(n.maxFreq);
        j.at("n_useOctaveDiv").get_to(n.useOctaveDiv);
        j.at("n_octaveDiv").get_to(n.octaveDiv);
        j.at("n_startSample").get_to(n.startSample);
        j.at("n_useList").get_to(n.useList);
    }

    // ---- OutInfo ----
    inline void from_json(const json& j, SihatFile::OutInfo& i) {
        j.at("outDir").get_to(i.outDir);
        j.at("prefix").get_to(i.prefix);
        j.at("extension").get_to(i.extension);
    }

    // ---- CorrelationSettings ----
    inline void from_json(const json& j, Sitrano::CorrelationSettings& c) {
        j.at("c_offset").get_to(c.periodStartOffsetMs);
        j.at("c_threshold").get_to(c.correlationThreshold);
    }

    // ---- OvertoneSettings ----
    inline void from_json(const json& j, Sitrano::OvertoneSettings& o) {
        j.at("o_useTolerance").get_to(o.useTolerance);
        j.at("o_tolerance").get_to(o.toleranceValue);
        j.at("o_postTransientStart").get_to(o.postTransientStart);
        j.at("o_FirstSample").get_to(o.overtoneFirstSample);
        j.at("o_useCustomSignal").get_to(o.useCustomSignal);
        j.at("o_sumAmplitudes").get_to(o.sumAmplitudes);
        j.at("o_mfft").get_to(o.fftSize);
        j.at("o_ignoreThreshold").get_to(o.overtoneThreshold);
        j.at("o_absThreshold").get_to(o.setAbsoluteThreshold);
    }

    // ---- HarmonicSettings ----
    inline void from_json(const json& j, Sitrano::HarmonicSettings& h) {
        j.at("h_enabled").get_to(h.applyHanning);
        j.at("h_Tolerance").get_to(h.toleranceValue);

        // Read the integer and cast it back to the enum type
        int styleAsInt;
        j.at("h_windowStyle").get_to(styleAsInt);
        h.style = static_cast<decltype(h.style)>(styleAsInt);
    }

    // ---- General Settings ----
    inline void from_json(const json& j, Sitrano::Settings& s) {
        j.at("fullHarmonicAnalysis").get_to(s.fullHarmonicAnalysis);
        j.at("pitchAnalysis").get_to(s.pitchAnalysis);
        j.at("transientSeparation").get_to(s.transientSeparation);
        j.at("periodAnalysis").get_to(s.periodAnalysis);
        j.at("overtoneAnalysis").get_to(s.overtoneAnalysis);
        j.at("harmonicAnalysis").get_to(s.harmonicAnalysis);
        j.at("noiseAnalysis").get_to(s.noiseAnalysis);
    }

    // ---- AnalysisConfig ----
    inline void from_json(const json& j, Sitrano::AnalysisConfig& c) {
        j.at("maxHarmonics").get_to(c.numHarmonics);
        j.at("N").get_to(c.nfft);
        j.at("hopSize").get_to(c.hopSize);
        j.at("startSample").get_to(c.startSample);
        j.at("tolerance").get_to(c.tolerance);
        j.at("verbose").get_to(c.verbose);
        j.at("bulkProcess").get_to(c.bulkProcess);

        from_json(j.at("pitchSettings"), c.pSettings);
        from_json(j.at("transientSettings"), c.tSettings);
        from_json(j.at("correlationSettings"), c.cSettings);
        from_json(j.at("overtoneSettings"), c.oSettings);
        from_json(j.at("harmonicSettings"), c.hSettings);
        from_json(j.at("noiseSettings"), c.nSettings);
    }


    inline void saveSettings(const Sitrano::AnalysisConfig& c,
        const Sitrano::Settings& s, const SihatFile::OutInfo& i, 
        const std::string& path) 
    { 
        json j; 
        json info, settings, config;

        to_json(info, i);
        to_json(settings, s);
        to_json(config, c);

        j["info"] = info; 
        j["settings"] = settings; 
        j["config"] = config; 
        std::ofstream file(path);
        if (!file.is_open()) 
        { 
            throw std::runtime_error("Could not open file for writing: " + path); 
        } 
        file << j.dump(4); 
        file.close(); 
        std::cout << "Configuration saved to: " << path << std::endl; 
    }

    inline void loadSettings(Sitrano::AnalysisConfig& c,
        Sitrano::Settings& s,
        SihatFile::OutInfo& i,
        const std::string& path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file for reading: " + path);
        }

        json j;
        try
        {
            file >> j; // Parse the JSON from the file
        }
        catch (json::parse_error& e)
        {
            file.close();
            throw std::runtime_error("Failed to parse JSON: " + std::string(e.what()));
        }
        file.close();

        // Check for the main keys
        if (!j.contains("info") || !j.contains("settings") || !j.contains("config"))
        {
            throw std::runtime_error("Invalid settings file format: " + path);
        }

        // Deserialize from the JSON object into your structs.
        // This will automatically use all the from_json functions we defined.
        try
        {
            from_json(j.at("info"), i);
            from_json(j.at("settings"), s);
            from_json(j.at("config"), c);
        }
        catch (json::exception& e)
        {
            throw std::runtime_error("Error reading settings from JSON: " + std::string(e.what()));
        }

        std::cout << "Configuration loaded from: " << path << std::endl;
    }
};