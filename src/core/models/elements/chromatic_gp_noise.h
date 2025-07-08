#pragma once

#include <unordered_set>
#include "../../../../eigen_config.h"
#include "../../types/model_element.h"

class chromatic_gp_noise_t : public model_element_t {
public:

    Eigen::VectorXd frequencies;
    int num_freqs;
    double days_per_coeff;
    double fixed_chromatic_idx;  // For Phase 1: fixed chromatic index

    chromatic_gp_noise_t(const std::optional<json_node_t>& element_json);
    
    // Calculate num_freqs based on time span
    void calculate_frequencies(double maxtspan);

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_valid_parameter(const string_t& param_name) const override;
    bool is_fully_specified() const override;
    string_t get_name() const override;

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    // Apply chromatic GP noise to power coefficient vector
    // Similar to red noise but with chromatic frequency scaling: (ν_ref/ν)^idx
    void apply(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const;

    // Check if idx is being fitted or fixed
    bool is_idx_fitted() const { 
        return get_optional_parameter("idx").has_value(); 
    }
    
    // Get the fixed idx value (when not fitted)
    double get_fixed_idx() const { 
        return fixed_chromatic_idx; 
    }

    // Calculate scaling correction factor for dynamic idx scaling
    // This accounts for the difference between current idx and reference idx used in matrix construction
    double calculate_scaling_correction(double current_idx, double reference_idx) const;

    void print() const override
    {
        std::cout << get_name() << " Element:" << std::endl;
        std::cout << "Number of Frequencies: " << num_freqs << std::endl;
        for (const auto& [name, param] : parameter_map_) {
            std::cout << name << ": ";
            param.print();
        }
    }

    void print_formatted(int element_index) const override
    {
        std::string extra_info = std::to_string(num_freqs) + " frequencies (" + std::to_string(days_per_coeff) + " days/coeff)";
        output_formatter::print_element_header(element_index, get_name(), extra_info);
        
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