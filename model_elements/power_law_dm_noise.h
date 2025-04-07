#pragma once

#include <unordered_set>
#include "../eigen_config.h"
#include "../types/model_element.h"

class pl_dm_noise_t : public model_element_t {
public:

    Eigen::VectorXd frequencies;
    int num_freqs;

    pl_dm_noise_t(const std::optional<json_node_t>& optional_config);

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
};