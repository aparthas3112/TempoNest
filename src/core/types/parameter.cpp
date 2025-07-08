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

    min_value = json_param.get_optional_value<double>("min_value").value_or(0);
    max_value = json_param.get_optional_value<double>("max_value").value_or(1);
}

void parameter_t::print() const
{
    std::cout << "Parameter: " << name << std::endl;
    std::cout << "Description: " << description << std::endl;
    std::cout << "Prior type: " << static_cast<int>(prior_type) << std::endl;
    std::cout << "Included: " << include << std::endl;
    std::cout << "Min value: " << min_value << std::endl;
    std::cout << "Max value: " << max_value << std::endl;
}

double parameter_t::get_value(const std::vector<double>& parameter_values) const
{
    if (index_ < 0 || index_ >= static_cast<int>(parameter_values.size())) {
        die("Invalid parameter index for parameter_t '" + name + "' (index: " + std::to_string(index_) + ")");
    }
    return parameter_values[index_];
}

double parameter_t::get_exp_value(const std::vector<double>& parameter_values) const
{
    if (index_ < 0 || index_ >= static_cast<int>(parameter_values.size())) {
        die("Invalid parameter index for parameter_t '" + name + "' (index: " + std::to_string(index_) + ")");
    }

    double exp_value = std::pow(10, parameter_values[index_]);
    return exp_value;
}