#pragma once

#include <iostream>
#include <vector>
#include "basic_types.h"
#include "parameter.h"

class model_element_t {
public:

    virtual ~model_element_t() = default;

    virtual bool is_valid_parameter(const string_t& param_name) const = 0;

    virtual bool is_fully_specified() const = 0;

    virtual string_t get_name() const = 0;

    virtual void set_parameter(const string_t& name, const parameter_t& param) = 0;

    virtual void print() const = 0;

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