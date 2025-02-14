#pragma once

#include "../eigen_config.h"
#include "../types/model_element.h"

class equad_t : public model_element_t {

private:

    Eigen::VectorXi flag_indices;

public:

    std::vector<string_t> flag_values;
    string_t flag;

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_fully_specified() const override;
    bool is_valid_parameter(const string_t& param_name) const override;
    string_t get_name() const override;

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    // function that applies the efac parameters to the noise vector
    void apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise, double& prior_term) const;
};