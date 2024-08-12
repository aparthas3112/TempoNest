#pragma once

#include <unordered_set>
#include <vector>
#include "../namespaces/settings.h"
#include "../types/model_element.h"
#include "../types/timing_parameter.h"

class timing_model_t : public model_element_t {
public:

    std::vector<timing_parameter_t> parameters;
    std::vector<string_t> t2_fitted_labels;
    std::vector<bool> marginalised;
    std::vector<long double> t2_fit_values;
    std::vector<long double> t2_fit_errors;
    std::vector<std::pair<int, int>> t2_fit_indices;

    int t2_total_fit;
    int design_size;

    // initialize all the vectors once the pulsar is loaded
    void initialise();

    // function to handle updating the residuals given the current set of parameters
    void update_residuals(double* Cube);

    void set_parameter(const string_t& name, const parameter_t& param) override;

    bool is_fully_specified() const override { return true; }

    bool is_valid_parameter(const string_t& param_name) const override { return true; }

    void print() const override
    {
        std::cout << "Timing Model Element:" << std::endl;
        for (int i = 0; i < parameters.size(); i++) {
            parameters[i].print();
        }
    }

    int get_fitted_dims() override { return t2_total_fit - design_size; }

    string_t get_name() const override { return "Timing Model"; }
};