#pragma once

#include <iostream>
#include <vector>
#include "basic_types.h"
#include "parameter.h"

// Forward declarations of all model element types
class model_element_t;
class efac_element;
class equad_element;
class power_law_red_noise_element;

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

#include "model_elements/efac.h"
#include "model_elements/equad.h"
#include "model_elements/power_law_red_noise.h"