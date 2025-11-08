#pragma once

#include"include/json.hpp"
#include"SihatJson.h"
#include"include/SitranoHeader.h"
#include"file/SihatFile.h"

using json = nlohmann::json;

namespace nlohmann {

    template <>
    struct adl_serializer<Sitrano::AnalysisConfig> {
        static void to_json(json& j, const Sitrano::AnalysisConfig& value) {
            SihatJson::to_json(j, value);
        }
    };

    template <>
    struct adl_serializer<Sitrano::PitchSettings> {
        static void to_json(json& j, const Sitrano::PitchSettings& value) {
            SihatJson::to_json(j, value);
        }
    };

    template<>
    struct adl_serializer<SihatFile::OutInfo> {
        static void to_json(json& j, const SihatFile::OutInfo& value) {
            SihatJson::to_json(j, value);
        }
    };
}