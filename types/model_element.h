#pragma once

#include <iostream>
#include <vector>
#include "basic_types.h"
#include "parameter.h"

class model_element_t {
public:

    virtual ~model_element_t() = default;
    virtual bool is_valid_parameter(const string_t& param_name) const = 0;
    virtual string_t get_name() const = 0;

    void add_parameter(const parameter_t& param)
    {
        if (is_valid_parameter(param.name)) {
            parameters.push_back(param);
        } else {
            std::cerr << "Warning: Ignoring invalid parameter '" << param.name << "' for "
                      << get_name() << std::endl;
        }
    }

    const std::vector<parameter_t>& get_parameters() const { return parameters; }

    void print() const
    {
        std::cout << get_name() << " parameters:" << std::endl;
        for (const auto& param : parameters) {
            param.print();
        }
    }

protected:

    std::vector<parameter_t> parameters;
};