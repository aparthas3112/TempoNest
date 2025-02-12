#pragma once

#include <memory>
#include <optional>
#include "../json/json_loader.h"
#include "../types/basic_types.h"
#include "../types/model_element.h"

class model_space_t {
public:

    explicit model_space_t();

    // Get element of specified type by name
    template <typename T>
    T& get_element(const string_t& name)
    {
        auto element = get_optional_element<T>(name);
        if (element) {
            return element->get();
        }
        die("Required element '" + name + "' not found in model");
    }

    /**
     * @brief Const version of get_element
     */
    template <typename T>
    const T& get_element(const string_t& name) const
    {
        auto element = get_optional_element<T>(name);
        if (element) {
            return element->get();
        }
        die("Required element '" + name + "' not found in model");
    }

    template <typename T>
    std::optional<std::reference_wrapper<T>> get_optional_element(const string_t& name)
    {
        auto it = elements_.find(name);
        if (it == elements_.end())
            return std::nullopt;

        auto* typed = it->second->as<T>();
        if (typed) {
            return std::ref(*typed);
        }
        die("Element '" + name + "' found with wrong type in model");
    }

    /**
     * @brief Const version of get_optional_element
     */
    template <typename T>
    std::optional<std::reference_wrapper<const T>> get_optional_element(const string_t& name) const
    {
        auto it = elements_.find(name);
        if (it == elements_.end())
            return std::nullopt;

        auto* typed = it->second->as<T>();
        if (typed) {
            return std::ref(*typed);
        }
        die("Element '" + name + "' found with wrong type in model");
    }

    // Get total number of fitted dimensions
    int get_fitted_dims() const;

    // Check if model is fully specified
    bool is_fully_specified() const;

    void update_array_size_info();

    void store_total_matrix();

    // Get matrices used for likelihood calculation
    const Eigen::MatrixXd& get_design_matrix() const { return design_matrix_; }

    const Eigen::MatrixXd& get_total_matrix() const { return total_matrix_; }

    // Get various sizes
    int get_total_size() const { return total_size_; }

    int get_design_size() const { return design_size_; }

    int get_noise_size() const { return noise_size_; }

    double get_max_tspan() const { return max_tspan_; }

    // Elements getter
    const std::unordered_map<string_t, element_t>& get_elements() const { return elements_; }

    /**
     * @brief Get parameters being sampled
     * @return Vector of pointers to non-fixed parameters that are included
     */
    std::vector<const parameter_t*> get_sampling_parameters() const;

private:

    // Helper methods
    void load_model();
    parameter_t parse_parameter(const json_node_t& json_param);
    element_t create_element(const string_t& type);

    void getEigenDVectorLike(Eigen::MatrixXd& design_matrix);

    // Member variables
    std::unordered_map<string_t, element_t> elements_;

    // Matrices and dimensions
    Eigen::MatrixXd design_matrix_;
    Eigen::MatrixXd total_matrix_;
    int total_size_{0};
    int design_size_{0};
    int noise_size_{0};
    double max_tspan_{0.0};
};