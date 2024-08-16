#pragma once

#include "../eigen_config.h"
#include "../types/model_element.h"

class efac_t : public model_element_t {
private:

    Eigen::VectorXi flag_indices;
    std::vector<string_t> flag_values;

public:

    optional_parameter_t global;
    optional_parameter_t per_flag;

    void set_parameter(const string_t& name, const parameter_t& param,
                       const rapidjson::Value& param_json) override;
    bool is_fully_specified() const override;
    bool is_valid_parameter(const string_t& param_name) const override;
    void print() const override;
    int get_fitted_dims() override;
    string_t get_name() const override;

    // function that applies the efac parameters to the noise vector
    void apply(double* Cube, Eigen::VectorXd& noise, double& prior_term, int& p_index);
};