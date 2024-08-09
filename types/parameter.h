#pragma once
#include "basic_types.h"

enum class prior_type_t { uniform, log_uniform, categorical };

class parameter_t {

public:

    string_t name;
    string_t description;
    prior_type_t prior_type;
    bool include;
    bool fit;
    double min_value;
    double max_value;

    void print() const
    {
        std::cout << "Parameter: " << name << std::endl;
        std::cout << "Description: " << description << std::endl;
        std::cout << "Prior type: " << static_cast<int>(prior_type) << std::endl;
        std::cout << "Included: " << include << std::endl;
        std::cout << "Fit: " << fit << std::endl;
        std::cout << "Min value: " << min_value << std::endl;
        std::cout << "Max value: " << max_value << std::endl;
    }
};