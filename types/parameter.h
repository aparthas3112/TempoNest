#pragma once
#include "basic_types.h"

enum class prior_type_t { uniform, log_uniform, categorical };

class parameter_t {

public:

    string_t name;
    string_t description;
    prior_type_t prior_type;
    bool is_included;
    bool should_fit;
    double min_value;
    double max_value;

    void print() const
    {
        std::cout << "Parameter: " << name << std::endl;
        std::cout << "Description: " << description << std::endl;
        std::cout << "Prior type: " << static_cast<int>(prior_type) << std::endl;
        std::cout << "Included: " << is_included << std::endl;
        std::cout << "Fit: " << should_fit << std::endl;
        std::cout << "Min value: " << min_value << std::endl;
        std::cout << "Max value: " << max_value << std::endl;
    }
};