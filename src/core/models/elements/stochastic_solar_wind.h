#pragma once

#include "../../../../eigen_config.h"
#include "../../types/model_element.h"

class stochastic_solar_wind_t : public model_element_t {
public:

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_fully_specified() const override;
    bool is_valid_parameter(const string_t& param_name) const override;
    string_t get_name() const override;

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    // Apply stochastic solar wind noise to noise vector
    // Formula: noise[i] += [10^(log_amp) × tdis2[i] / ne_sw]²
    void apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise) const;

    void print_formatted(int element_index) const override
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
};