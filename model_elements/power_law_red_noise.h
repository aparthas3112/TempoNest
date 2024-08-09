#pragma once

#include <unordered_set>
#include "../types/model_element.h"

class red_noise_element : public model_element_t {
public:

    bool is_valid_parameter(const string_t& param_name) const override
    {
        static const std::unordered_set<string_t> valid_params = {"amplitude", "spectral_index"};
        return valid_params.find(param_name) != valid_params.end();
    }

    string_t get_name() const override { return "Power Law Red Noise"; }
};