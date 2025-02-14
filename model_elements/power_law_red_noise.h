#pragma once

#include <unordered_set>
#include "../eigen_config.h"
#include "../types/model_element.h"

class pl_red_noise_t : public model_element_t {
public:

    Eigen::VectorXd frequencies;
    int num_freqs;

    pl_red_noise_t();

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_valid_parameter(const string_t& param_name) const override;
    bool is_fully_specified() const override;
    string_t get_name() const override;

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    // use the red noise params to calculate the red noise power at each frequency
    void apply(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const;
};