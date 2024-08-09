#pragma once

#include <memory>
#include "rapidjson/document.h"
#include "types/model_element.h"
#include "types/parameter.h"

class json_loader {
public:

    static parameter_t parse_parameter(const json_value_t& json_param);
    static element_t create_model_element(const string_t& element_name);
    static void load_from_json(const string_t& filename);

private:

    static rapidjson::Document parse_json_file(const string_t& filename);
};