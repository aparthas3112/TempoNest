#pragma once
#include <string>
#include "../json/json_loader.h"

enum class prior_type_t { uniform, log_uniform, categorical };

class parameter_t {
public:

    std::string name;
    std::string description;
    prior_type_t prior_type;
    bool include;
    bool fit;
    double min_value;
    double max_value;

    parameter_t();
    virtual ~parameter_t() = default;

    virtual void load_from_json(const json_node_t& json_param);
    virtual void print() const;
    virtual double get_value(double& cube_value) const;
    virtual double get_exp_value(double& cube_value) const;
};