#pragma once

#include <memory>
#include <optional>
#include "../json/json_loader.h"
#include "../model/model_space.h"
#include "basic_types.h"
#include "model_element.h"

namespace model {

extern model_space_t model_space;

void load_model(const string_t& filename);

// Wrapper functions that delegate to model_space (optional, for backwards compatibility)
template <typename T>
T& get_element(const string_t& name)
{
    return model_space.get_element<T>(name);
}

template <typename T>
std::optional<std::reference_wrapper<T>> get_optional_element(const string_t& name)
{
    return model_space.get_optional_element<T>(name);
}

}  // namespace model