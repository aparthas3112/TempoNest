#pragma once

#include <unordered_set>
#include "../types/model_element.h"

class equad_element : public model_element_t {
public:

    optional_parameter_t global;

    void set_parameter(const string_t& name, const parameter_t& param) override
    {
        if (name == "global") {
            global = param;
        } else {
            throw std::runtime_error("Invalid parameter name for EQUAD: " + name);
        }
    }

    bool is_fully_specified() const override { return true; }

    bool is_valid_parameter(const string_t& param_name) const override
    {
        static const std::unordered_set<string_t> valid_params = {"global"};
        return valid_params.find(param_name) != valid_params.end();
    }

    void print() const override
    {
        std::cout << "EQUAD Element:" << std::endl;
        if (global.has_value()) {
            std::cout << "global: ";
            global->print();
        }
    }

    int get_fitted_dims() override
    {
        if (global.has_value()) {
            return 1;
        }
        return 0;
    }

    string_t get_name() const override { return "EQUAD"; }
};