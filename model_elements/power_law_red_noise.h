#pragma once

#include <unordered_set>
#include "../eigen_config.h"
#include "../types/model_element.h"

class pl_red_noise_element : public model_element_t {
public:

    parameter_t amplitude;
    parameter_t spectral_index;
    optional_parameter_t num_coeffs;
    Eigen::VectorXd frequencies;

    pl_red_noise_element()
    {
        frequencies = Eigen::VectorXd::Zero(33);
        for (int i = 0; i < frequencies.size(); i++) {
            frequencies[i] = i + 1;
        }
    }

    void set_parameter(const string_t& name, const parameter_t& param) override
    {
        if (name == "amplitude") {
            amplitude = param;
        } else if (name == "spectral_index") {
            spectral_index = param;
        } else if (name == "num_coeffs") {
            num_coeffs = param;
        } else {
            throw std::runtime_error("Invalid parameter name for Power Law Red Noise: " + name);
        }
    }

    bool is_valid_parameter(const string_t& param_name) const override
    {
        static const std::unordered_set<string_t> valid_params = {"amplitude", "spectral_index",
                                                                  "num_coeffs"};
        return valid_params.find(param_name) != valid_params.end();
    }

    bool is_fully_specified() const override
    {
        // If either amplitude or spectral_index is included, both must be
        if (amplitude.include == spectral_index.include) {
            return true;
        }
        if (num_coeffs.has_value() && num_coeffs.value().include && amplitude.include) {
            return true;
        }

        return false;
    }

    void print() const override
    {
        std::cout << "Power Law Red Noise Element:" << std::endl;
        std::cout << "Amplitude: ";
        amplitude.print();
        std::cout << "Spectral Index: ";
        spectral_index.print();
        if (num_coeffs.has_value()) {
            std::cout << "Cutoff: ";
            num_coeffs->print();
        }
    }

    string_t get_name() const override { return "Power Law Red Noise"; }
};