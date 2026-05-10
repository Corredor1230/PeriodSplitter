#pragma once

#include"include/json.hpp"
#include"include/SitranoHeader.h"
#include"file/SihatFile.h"
#include<fstream>

using json = nlohmann::json;

namespace SihatJson {
    // ---- PitchSettings ----
    inline void to_json(json& j, const Sihat::PitchSettings& p) {
        j = json{
            {"p_threshold", p.modeThreshold},
            {"p_inCents", p.toleranceInCents},
            {"p_minFreq", p.minFreq},
            {"p_maxFreq", p.maxFreq}
        };
    }

    // ---- TransientSettings ----
    inline void to_json(json& j, const Sihat::TransientSettings& t) {
        j = json{
            {"t_StartSample", t.tStartSample},
            {"t_UseMs", t.useMs},
            {"t_rmsSampleSize", t.rmsSampleSize},
            {"t_rmsHopLength", t.rmsSampleHopLength},
            {"t_rmsSizeMs", t.transientRmsSizeMs},
            {"t_rmsHopRatio", t.transientRmsHopRatio},
            {"t_factor", t.transientFactor},
            {"t_threshold", t.transientThreshold},
            {"t_preAttack", t.preAttack},
            {"t_correlationThreshold", t.tCorrelationThreshold}
        };
    }

    //Single transient settings
    inline void to_json(json& j, const Sihat::SingleTransientSettings& st){
        j = json{
            {"st_rmsWindow", st.rmsWindow},
            {"st_rmsHopSize", st.rmsHopSize},
            {"st_nfft", st.nfft},
            {"st_hopSize", st.hopSize},
            {"st_inThreshold", st.inThreshold},
            {"st_outThreshold", st.outThreshold},
            {"st_numBands", st.numBands},
            {"st_maxOvertones", st.maxOvertones},
            {"st_startSample", st.startSample}
        };
    }

    inline void to_json(json& j, const Sihat::TransientFFTSettings& tfft){
        j = json{
            {"tfft_nfft", tfft.nfft},
            {"tfft_hopSize", tfft.hopSize},
            {"tfft_useFFT", tfft.useFFT},
            {"tfft_flatnessThreshold", tfft.flatnessThreshold}
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
    inline void to_json(json& j, const Sihat::CorrelationSettings& c) {
        j = json{
            {"c_offset", c.periodStartOffsetMs},
            {"c_threshold", c.correlationThreshold}
        };
    }

    // ---- OvertoneSettings ----
    inline void to_json(json& j, const Sihat::OvertoneSettings& o) {
        j = json{
            {"o_useTolerance", o.useTolerance},
            {"o_useTolInHz", o.useTolInHz},
            {"o_tolInHz", o.tolInHz},
            {"o_chooseFirstSample", o.chooseFirstSample},
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
    inline void to_json(json& j, const Sihat::HarmonicSettings& h) {
        j = json{
            {"h_enabled", h.applyHanning},
            {"h_windowStyle", static_cast<int>(h.style)},  // convert enum to int
            {"h_Tolerance", h.toleranceValue},
            {"h_tolInHz", h.tolInHz},
            {"h_useTolInHz", h.useTolInHz}
        };
    }

    // --- NoiseSettings ---
    inline void to_json(json& j, const Sihat::NoiseSettings& n) {
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

    // ---- HPSSSettings ----
    inline void to_json(json& j, const Sihat::HPSSSettings& hp) {
        j = json{
            {"hpss_nfft", hp.nfft},
            {"hpss_hopSize", hp.hopSize},
            {"hpss_filtSize", hp.filtSize}
        };
    }

    // ---- General Settings ----
    inline void to_json(json& j, const Sihat::Settings& s) {
        j = json{
            {"fullHarmonicAnalysis", s.fullHarmonicAnalysis},
            {"pitchAnalysis", s.pitchAnalysis},
            {"sourceSeparation", s.sourceSeparation},
            {"transientSeparation", s.transientSeparation},
            {"periodAnalysis", s.periodAnalysis},
            {"overtoneAnalysis", s.overtoneAnalysis},
            {"harmonicAnalysis", s.harmonicAnalysis},
            {"noiseAnalysis", s.noiseAnalysis}
        };
    }

    // ---- AnalysisConfig ----
    inline void to_json(json& j, const Sihat::AnalysisConfig& c) {

        json pSettings;
        json tSettings;
        json hpssSettings;
        json tfftSettings;
        json stSettings;
        json cSettings;
        json oSettings;
        json hSettings;
        json nSettings;

        to_json(pSettings, c.pSettings);
        to_json(tSettings, c.tSettings);
        to_json(tfftSettings, c.tfftSettings);
        to_json(cSettings, c.cSettings);
        to_json(oSettings, c.oSettings);
        to_json(hSettings, c.hSettings);
        to_json(nSettings, c.nSettings);
        to_json(hpssSettings, c.hpSettings);
        to_json(stSettings, c.stSettings);


        j = json{
            {"maxHarmonics", c.numHarmonics},
            {"N", c.nfft},
            {"hopSize", c.hopSize},
            {"startSample", c.startSample},
            {"tolerance", c.tolerance},
            {"hpssSettings", hpssSettings},
            {"pitchSettings", pSettings},
            {"transientSettings", tSettings},
            {"stSettings", stSettings},
            {"transientFFTSettings", tfftSettings},
            {"correlationSettings", cSettings},
            {"overtoneSettings", oSettings},
            {"harmonicSettings", hSettings},
            {"noiseSettings", nSettings},
            {"verbose", c.verbose},
            {"bulkProcess", c.bulkProcess}
        };
    }

    // --- PitchSettings ---
    inline void from_json(const json& j, Sihat::PitchSettings& p) {
        j.at("p_threshold").get_to(p.modeThreshold);
        j.at("p_inCents").get_to(p.toleranceInCents);
        j.at("p_minFreq").get_to(p.minFreq);
        j.at("p_maxFreq").get_to(p.maxFreq);
    }

    // ---- TransientSettings ----
    inline void from_json(const json& j, Sihat::TransientSettings& t) {
        j.at("t_StartSample").get_to(t.tStartSample);
        j.at("t_UseMs").get_to(t.useMs);
        j.at("t_rmsSampleSize").get_to(t.rmsSampleSize);
        j.at("t_rmsHopLength").get_to(t.rmsSampleHopLength);
        j.at("t_rmsSizeMs").get_to(t.transientRmsSizeMs);
        j.at("t_rmsHopRatio").get_to(t.transientRmsHopRatio);
        j.at("t_factor").get_to(t.transientFactor);
        j.at("t_threshold").get_to(t.transientThreshold);
        j.at("t_preAttack").get_to(t.preAttack);
        j.at("t_correlationThreshold").get_to(t.tCorrelationThreshold);
    }

    inline void from_json(const json& j, Sihat::TransientFFTSettings& tfft){
        j.at("tfft_nfft").get_to(tfft.nfft);
        j.at("tfft_hopSize").get_to(tfft.hopSize);
        j.at("tfft_useFFT").get_to(tfft.useFFT);
        j.at("tfft_flatnessThreshold").get_to(tfft.flatnessThreshold);
    }

    //Single transient settings
    inline void from_json(const json& j, Sihat::SingleTransientSettings& st)
    {
        j.at("st_rmsWindow").get_to(st.rmsWindow);
        j.at("st_rmsHopSize").get_to(st.rmsHopSize); 
        j.at("st_nfft").get_to(st.nfft);
        j.at("st_hopSize").get_to(st.hopSize);
        j.at("st_inThreshold").get_to(st.inThreshold);
        j.at("st_outThreshold").get_to(st.outThreshold);
        j.at("st_numBands").get_to(st.numBands);
        j.at("st_startSample").get_to(st.startSample);
        j.at("st_maxOvertones").get_to(st.maxOvertones);
    }

    // --- NoiseSettings ---
    inline void from_json(const json& j, Sihat::NoiseSettings& n) {
        j.at("n_nfft").get_to(n.nfft);
        j.at("n_hopSize").get_to(n.hopSize);
        j.at("n_numBins").get_to(n.numBins);
        j.at("n_minFreq").get_to(n.minFreq);
        j.at("n_maxFreq").get_to(n.maxFreq);
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
    inline void from_json(const json& j, Sihat::CorrelationSettings& c) {
        j.at("c_offset").get_to(c.periodStartOffsetMs);
        j.at("c_threshold").get_to(c.correlationThreshold);
    }

    // ---- OvertoneSettings ----
    inline void from_json(const json& j, Sihat::OvertoneSettings& o) {
        j.at("o_useTolerance").get_to(o.useTolerance);
        j.at("o_useTolInHz").get_to(o.useTolInHz);
        j.at("o_tolInHz").get_to(o.tolInHz);
        j.at("o_chooseFirstSample").get_to(o.chooseFirstSample);
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
    inline void from_json(const json& j, Sihat::HarmonicSettings& h) {
        j.at("h_enabled").get_to(h.applyHanning);
        j.at("h_Tolerance").get_to(h.toleranceValue);
        j.at("h_tolInHz").get_to(h.tolInHz);
        j.at("h_useTolInHz").get_to(h.useTolInHz);

        // Read the integer and cast it back to the enum type
        int styleAsInt;
        j.at("h_windowStyle").get_to(styleAsInt);
        h.style = static_cast<decltype(h.style)>(styleAsInt);
    }

    // ---- HPSSSettings ----
    inline void from_json(const json& j, Sihat::HPSSSettings& hp) {
        j.at("hpss_nfft").get_to(hp.nfft);
        j.at("hpss_hopSize").get_to(hp.hopSize);
        j.at("hpss_filtSize").get_to(hp.filtSize);
    }

    // ---- General Settings ----
    inline void from_json(const json& j, Sihat::Settings& s) {
        j.at("fullHarmonicAnalysis").get_to(s.fullHarmonicAnalysis);
        j.at("sourceSeparation").get_to(s.sourceSeparation);
        j.at("pitchAnalysis").get_to(s.pitchAnalysis);
        j.at("transientSeparation").get_to(s.transientSeparation);
        j.at("periodAnalysis").get_to(s.periodAnalysis);
        j.at("overtoneAnalysis").get_to(s.overtoneAnalysis);
        j.at("harmonicAnalysis").get_to(s.harmonicAnalysis);
        j.at("noiseAnalysis").get_to(s.noiseAnalysis);
    }

    // ---- AnalysisConfig ----
    inline void from_json(const json& j, Sihat::AnalysisConfig& c) {
        j.at("maxHarmonics").get_to(c.numHarmonics);
        j.at("N").get_to(c.nfft);
        j.at("hopSize").get_to(c.hopSize);
        j.at("startSample").get_to(c.startSample);
        j.at("tolerance").get_to(c.tolerance);
        j.at("verbose").get_to(c.verbose);
        j.at("bulkProcess").get_to(c.bulkProcess);

        from_json(j.at("pitchSettings"), c.pSettings);
        from_json(j.at("transientSettings"), c.tSettings);
        from_json(j.at("hpssSettings"), c.hpSettings);
        from_json(j.at("transientFFTSettings"), c.tfftSettings);
        from_json(j.at("stSettings"), c.stSettings);
        from_json(j.at("correlationSettings"), c.cSettings);
        from_json(j.at("overtoneSettings"), c.oSettings);
        from_json(j.at("harmonicSettings"), c.hSettings);
        from_json(j.at("noiseSettings"), c.nSettings);
    }


    inline void saveSettings(const Sihat::AnalysisConfig& c,
        const Sihat::Settings& s, const SihatFile::OutInfo& i, 
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

    inline void loadSettings(Sihat::AnalysisConfig& c,
        Sihat::Settings& s,
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