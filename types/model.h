#pragma once

#include <memory>
#include <optional>
#include "basic_types.h"
#include "model_element.h"

namespace model {

// model elements
extern optional_element_t pl_red_noise;
extern optional_element_t pl_dm_noise;
extern optional_element_t efac;
extern optional_element_t equad;

void load_model(const string_t& filename);
}  // namespace model