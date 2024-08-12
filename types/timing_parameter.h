#pragma once
#include <math.h> /* pow */
#include <iostream>
#include "basic_types.h"
#include "parameter.h"

class timing_parameter_t : public parameter_t {
public:

    long double ld_pmin = 0;
    long double ld_pmax = 0;
    int t2_p_index;
    int t2_k_index;

    // Constructor that takes a parameter_t
    timing_parameter_t(const parameter_t& param) : parameter_t(param) {}

    void print() const override
    {
        parameter_t::print();
        std::cout << "LD Min value: " << ld_pmin << std::endl;
        std::cout << "LD Max value: " << ld_pmax << std::endl;
    }
};