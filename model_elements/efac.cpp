#include "efac.h"
#include <iostream>

void efac_t::set_parameter(const string_t& name, const parameter_t& param)
{
    if (name == "global") {
        global = param;
    } else {
        throw std::runtime_error("Invalid parameter name for EFAC: " + name);
    }
}

bool efac_t::is_fully_specified() const
{
    return true;
}

bool efac_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"global"};
    return valid_params.find(param_name) != valid_params.end();
}

void efac_t::print() const
{
    std::cout << "EFAC Element:" << std::endl;
    if (global.has_value()) {
        std::cout << "global: ";
        global->print();
    }
}

int efac_t::get_fitted_dims()
{
    if (global.has_value()) {
        return 1;
    }
    return 0;
}

string_t efac_t::get_name() const
{
    return "EFAC";
}

void efac_t::apply(double* Cube, Eigen::VectorXd& noise, double& prior_term, int& p_index)
{
    if (global.has_value()) {
        double multiplier = global.value().get_exp_value(Cube[p_index]);

        if (global.value().prior_type == prior_type_t::uniform) {
            prior_term += log(multiplier);
        }

        p_index++;

        noise = (noise * multiplier).array().square();
    }
}