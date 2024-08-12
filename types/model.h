#pragma once

#include <memory>
#include <optional>
#include "basic_types.h"
#include "model_element.h"

namespace model {

extern pulsar_t* pulsar;

extern double max_tspan;

extern int total_size;
extern int design_size;
extern int noise_size;
extern Eigen::MatrixXd design_matrix;
extern Eigen::MatrixXd total_matrix;

// model elements
extern optional_element_t pl_red_noise;
extern optional_element_t pl_dm_noise;
extern optional_element_t efac;
extern optional_element_t equad;

void load_model(const string_t& filename);
}  // namespace model