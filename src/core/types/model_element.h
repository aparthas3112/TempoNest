#pragma once

#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <vector>
#include "basic_types.h"
#include "parameter.h"
#include "rapidjson/document.h"
#include "../utils/output_formatter.h"

class model_element_t {
protected:

    std::unordered_map<string_t, parameter_t> parameter_map_;

public:

    /**
     * @brief Get read-only access to the parameters vector
     * @return Const reference to parameters
     */
    std::vector<const parameter_t*> get_parameters() const
    {
        std::vector<const parameter_t*> params;
        params.reserve(parameter_map_.size());
        for (const auto& [key, param] : parameter_map_) {
            const_cast<parameter_t&>(param).set_id(key);  // Need const_cast since we're in a const function
            params.push_back(&param);
        }
        return params;
    }

    virtual ~model_element_t() = default;

    virtual bool is_valid_parameter(const string_t& param_name) const = 0;

    virtual bool is_fully_specified() const = 0;

    virtual string_t get_name() const = 0;

    virtual void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) = 0;

    int get_fitted_dims() const { return parameter_map_.size(); }

    virtual void print() const
    {
        std::cout << get_name() << " Element:" << std::endl;
        for (const auto& [name, param] : parameter_map_) {
            std::cout << name << ": ";
            param.print();
        }
    }

    virtual void print_formatted(int element_index) const
    {
        output_formatter::print_element_header(element_index, get_name());
        
        if (parameter_map_.empty()) {
            output_formatter::print_marginalised_note();
            return;
        }
        
        size_t param_count = 0;
        size_t total_params = parameter_map_.size();
        
        for (const auto& [name, param] : parameter_map_) {
            param_count++;
            bool is_last = (param_count == total_params);
            
            std::string prior_str = (param.prior_type == prior_type_t::uniform) ? "uniform" : "log-uniform";
            output_formatter::print_parameter(name, prior_str, param.min_value, param.max_value, param.include, is_last);
        }
    }

    virtual void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const = 0;

    template <typename T>
    T* as()
    {
        return dynamic_cast<T*>(this);
    }

    template <typename T>
    const T* as() const
    {
        return dynamic_cast<const T*>(this);
    }

    std::optional<const parameter_t*> get_optional_parameter(const string_t& name) const
    {
        auto it = parameter_map_.find(name);
        if (it == parameter_map_.end()) {
            return std::nullopt;
        }
        return &it->second;
    }

    const parameter_t* get_parameter(const string_t& name) const
    {
        auto param = get_optional_parameter(name);
        if (!param) {
            die("Parameter not found: " + name);
        }
        return *param;
    }

    const std::unordered_map<string_t, parameter_t>& get_parameter_map() const { return parameter_map_; }
};

// we need to include all the model elements here to avoid circular dependencies
#include "../models/elements/efac.h"
#include "../models/elements/equad.h"
#include "../models/elements/ecorr.h"
#include "../models/elements/power_law_dm_noise.h"
#include "../models/elements/power_law_red_noise.h"
#include "../models/elements/timing_model.h"
#include "../models/elements/deterministic_solar_wind.h"
#include "../models/elements/stochastic_solar_wind.h"
#include "../models/elements/chromatic_gp_noise.h"