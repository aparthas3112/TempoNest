#include "model.h"
#include "../json_loader.h"

namespace model {

pulsar_t* pulsar;

int total_size;
int design_size;
int noise_size;
double max_tspan;

// matrices, can be used to speed up likelihood if they are constant
Eigen::MatrixXd design_matrix;
Eigen::MatrixXd total_matrix;

optional_element_t pl_red_noise;
optional_element_t pl_dm_noise;
optional_element_t efac;
optional_element_t equad;

void load_model(const string_t& filename)
{
    json_loader::load_model(filename);
}

}  // namespace model