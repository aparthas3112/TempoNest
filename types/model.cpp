#include "model.h"
#include "../json_loader.h"

namespace model {

int total_size;
int design_size;
int noise_size;
double max_tspan;

// matrices, can be used to speed up likelihood if they are constant
Eigen::MatrixXd design_matrix;
Eigen::MatrixXd total_matrix;

element_t timing_model;

optional_element_t pl_red_noise;
optional_element_t pl_dm_noise;
optional_element_t efac;
optional_element_t equad;

void load_model(const string_t& filename)
{
    json_loader::load_model(filename);
}

int get_model_dims()
{
    int dims = 0;

    if (pl_red_noise.has_value()) {
        dims += pl_red_noise.value()->get_fitted_dims();
    }

    if (pl_dm_noise.has_value()) {
        dims += pl_dm_noise.value()->get_fitted_dims();
    }

    if (efac.has_value()) {
        dims += efac.value()->get_fitted_dims();
    }

    if (equad.has_value()) {
        dims += equad.value()->get_fitted_dims();
    }

    dims += timing_model->get_fitted_dims();

    return dims;
}
}  // namespace model