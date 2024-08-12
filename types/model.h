#pragma once

#include <memory>
#include <optional>
#include "basic_types.h"
#include "model_element.h"

namespace model {

extern double max_tspan;

extern int total_size;
extern int design_size;
extern int noise_size;
extern Eigen::MatrixXd design_matrix;
extern Eigen::MatrixXd total_matrix;

// model elements

// there is always some kind of timing model even if we are just marginalising over it
extern element_t timing_model;

// all other properties are optional
extern optional_element_t pl_red_noise;
extern optional_element_t pl_dm_noise;
extern optional_element_t efac;
extern optional_element_t equad;

void load_model(const string_t& filename);

int get_model_dims();
}  // namespace model