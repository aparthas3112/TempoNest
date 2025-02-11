#include "parameter.h"
#include <cmath>
#include <iostream>

parameter_t::parameter_t() : min_value(0), max_value(1) {}

void parameter_t::load_from_json(const json_node_t& json_param)
{
    name = json_param.get_value<string_t>("name");
    description = json_param.get_optional_value<string_t>("description").value_or("");

    prior_type = string_t(json_param.get_value<string_t>("prior_type")) == "uniform" ? prior_type_t::uniform : prior_type_t::log_uniform;

    include = json_param.get_value<bool>("include");
    fit = json_param.get_value<bool>("fit");

    min_value = json_param.get_optional_value<double>("min_value").value_or(0);
    max_value = json_param.get_optional_value<double>("max_value").value_or(1);
}

void parameter_t::print() const
{
    std::cout << "Parameter: " << name << std::endl;
    std::cout << "Description: " << description << std::endl;
    std::cout << "Prior type: " << static_cast<int>(prior_type) << std::endl;
    std::cout << "Included: " << include << std::endl;
    std::cout << "Fit: " << fit << std::endl;
    std::cout << "Min value: " << min_value << std::endl;
    std::cout << "Max value: " << max_value << std::endl;
}

double parameter_t::get_value(double& cube_value) const
{
    double scaled_value = (max_value - min_value) * cube_value + min_value;
    cube_value = scaled_value;
    return scaled_value;
}

double parameter_t::get_exp_value(double& cube_value) const
{
    double scaled_value = (max_value - min_value) * cube_value + min_value;
    cube_value = scaled_value;
    scaled_value = std::pow(10, scaled_value);
    return scaled_value;
}