#pragma once

#include <unordered_set>
#include "../eigen_config.h"
#include "../types/model_element.h"

class pl_dm_noise_t : public model_element_t {
public:

    parameter_t amplitude;
    parameter_t spectral_index;
    optional_parameter_t num_coeffs;
    Eigen::VectorXd frequencies;
    int num_freqs;

    pl_dm_noise_t();

    void set_parameter(const string_t& name, const parameter_t& param,
                       const rapidjson::Value& param_json) override;
    bool is_valid_parameter(const string_t& param_name) const override;
    bool is_fully_specified() const override;
    void print() const override;
    int get_fitted_dims() override;
    string_t get_name() const override;

    // use the dm noise params to calculate the dm noise power at each frequency
    void apply(double* Cube, Eigen::VectorXd& powercoeff, int& p_count, int& start_pos,
               double maxtspan, double& uniform_prior, double& freq_det);
};