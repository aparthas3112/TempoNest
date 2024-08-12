#include "parameter.h"
#include <cmath>
#include <iostream>
#include "../json_loader.h"

parameter_t::parameter_t() : min_value(0), max_value(1) {}

void parameter_t::load_from_json(const rapidjson::Value& json_param)
{
    name = json_param["name"].GetString();
    json_loader::get_if_present(json_param, "description", description);

    prior_type = string_t(json_param["prior_type"].GetString()) == "uniform"
                     ? prior_type_t::uniform
                     : prior_type_t::log_uniform;

    include = json_param["include"].GetBool();
    fit = json_param["fit"].GetBool();

    // default priors are just [0,1]
    json_loader::get_if_present(json_param, "min_value", min_value);
    json_loader::get_if_present(json_param, "max_value", max_value);
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