#pragma once
#include <optional>
#include <string>
#include "rapidjson/document.h"
#include "tempo2.h"

// Forward declaration
class model_element_t;
class parameter_t;
class pulsar;

using pulsar_t = pulsar;

using string_t = std::string;
using json_value_t = rapidjson::Value;

using element_t = std::unique_ptr<model_element_t>;
using optional_element_t = std::optional<element_t>;

using optional_parameter_t = std::optional<parameter_t>;
