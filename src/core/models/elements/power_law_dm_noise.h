#pragma once

#include <unordered_set>
#include "../../../../eigen_config.h"
#include "../../types/model_element.h"

class pl_dm_noise_t : public model_element_t {
public:

    Eigen::VectorXd frequencies;
    int num_freqs;
    double days_per_coeff;

    pl_dm_noise_t(const std::optional<json_node_t>& optional_config);
    
    // Calculate num_freqs based on time span
    void calculate_frequencies(double maxtspan);

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_valid_parameter(const string_t& param_name) const override;
    bool is_fully_specified() const override;
    string_t get_name() const override;

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    // use the dm noise params to calculate the dm noise power at each frequency
    void apply(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const;

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