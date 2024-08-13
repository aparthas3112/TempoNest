#include "equad.h"
#include <iostream>

void equad_t::set_parameter(const string_t& name, const parameter_t& param)
{
    if (name == "global") {
        global = param;
    } else {
        throw std::runtime_error("Invalid parameter name for EQUAD: " + name);
    }
}

bool equad_t::is_fully_specified() const
{
    return true;
}

bool equad_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"global"};
    return valid_params.find(param_name) != valid_params.end();
}

void equad_t::print() const
{
    std::cout << "EQUAD Element:" << std::endl;
    if (global.has_value()) {
        std::cout << "global: ";
        global->print();
    }
}

int equad_t::get_fitted_dims()
{
    if (global.has_value()) {
        return 1;
    }
    return 0;
}

string_t equad_t::get_name() const
{
    return "EQUAD";
}

void equad_t::apply(double* Cube, Eigen::VectorXd& noise, double& prior_term, int& p_index)
{
    if (global.has_value()) {
        double value = global.value().get_exp_value(Cube[p_index]);
        double equad = value * value;

        if (global.value().prior_type == prior_type_t::uniform) {
            prior_term += 0.5 * log(value);
        }

        p_index++;

        noise = noise.array() + equad;
    }
}