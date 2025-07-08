#pragma once
#include <string>
#include "../../json/json_loader.h"

enum class prior_type_t { uniform, log_uniform, categorical };

class parameter_t {
public:

    std::string name;
    std::string id = "";
    std::string description;
    prior_type_t prior_type;
    bool include;
    double min_value;
    double max_value;

    parameter_t();
    virtual ~parameter_t() = default;

    virtual void load_from_json(const json_node_t& json_param);
    virtual void print() const;
    virtual double get_value(const std::vector<double>& parameter_values) const;
    virtual double get_exp_value(const std::vector<double>& parameter_values) const;

    /**
     * @brief Get parameter index in sampler's parameter vector
     *
     * @return int Index (-1 if parameter is not being sampled)
     */
    int get_index() const { return index_; }

    /**
     * @brief Set parameter index in sampler's parameter vector
     *
     * @param index New index value
     */
    void set_index(int index) { index_ = index; }

    void set_parent(model_element_t* parent) { parent_ = parent; }

    model_element_t* get_parent() const { return parent_; }

    void set_id(std::string id) { this->id = id; }

private:

    model_element_t* parent_{nullptr};
    /**
     * Index in sampler's parameter vector.
     * -1 indicates parameter is not being sampled.
     */
    int index_ = -1;
};