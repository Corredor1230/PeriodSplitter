#pragma once

#include"include/json.hpp"
#include"SihatJson.h"
#include"include/SitranoHeader.h"
#include"file/SihatFile.h"

using json = nlohmann::json;

namespace nlohmann {

    template <>
    struct adl_serializer<Sihat::AnalysisConfig> {
        static void to_json(json& j, const Sihat::AnalysisConfig& value) {
            SihatJson::to_json(j, value);
        }
    };

    template <>
    struct adl_serializer<Sihat::PitchSettings> {
        static void to_json(json& j, const Sihat::PitchSettings& value) {
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