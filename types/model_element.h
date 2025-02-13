#pragma once

#include <cstdio>
#include <iostream>
#include <vector>
#include "basic_types.h"
#include "parameter.h"
#include "rapidjson/document.h"

class model_element_t {
protected:

    /// Vector of pointers to parameters used by this element
    std::vector<const parameter_t*> parameters_;

public:

    /**
     * @brief Get read-only access to the parameters vector
     * @return Const reference to parameters
     */
    const std::vector<const parameter_t*>& get_parameters() const { return parameters_; }

    virtual ~model_element_t() = default;

    virtual bool is_valid_parameter(const string_t& param_name) const = 0;

    virtual bool is_fully_specified() const = 0;

    virtual string_t get_name() const = 0;

    virtual void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) = 0;

    virtual int get_fitted_dims() = 0;

    virtual void print() const = 0;

    virtual void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const = 0;

    template <typename T>
    T* as()
    {
        return dynamic_cast<T*>(this);
    }

    template <typename T>
    const T* as() const
    {
        return dynamic_cast<const T*>(this);
    }
};

// we need to include all the model elements here to avoid circular dependencies
#include "model_elements/efac.h"
#include "model_elements/equad.h"
#include "model_elements/power_law_dm_noise.h"
#include "model_elements/power_law_red_noise.h"
#include "model_elements/timing_model.h"