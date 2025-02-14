#include <TempoNest.h>
#include "../logger.h"
#include "tests.h"

bool run_likelihood_tests(const std::shared_ptr<model_t> model)
{

    std::cout << "get the model space and likelihood" << std::endl;
    const model_space_t& model_space = model->get_model_space();
    const likelihood_t& likelihood = model->get_likelihood();

    std::cout << "get dims" << std::endl;

    int ndims = model_space.get_fitted_dims();
    int npars = ndims;

    double* Cube = new double[5];
    Cube[0] = 0.91;
    Cube[1] = 0.85;
    Cube[2] = 0.85;
    Cube[3] = 0.9;
    Cube[4] = 0.8;

    std::cout << "get params" << std::endl;

    std::vector<const parameter_t*> parameters = model->get_sampling_parameters();

    std::cout << "Have params " << parameters.size() << " " << ndims << std::endl;
    std::vector<double> params(ndims);

    for (int i = 0; i < ndims; ++i) {
        double physical_value = parameters[i]->min_value + Cube[i] * (parameters[i]->max_value - parameters[i]->min_value);
        params[i] = physical_value;
    }

    std::cout << "run likelihood" << std::endl;

    double result = likelihood(model_space, params);

    if (fabs(result - 956.0983174) > 1e-6) {
        logger::log_error("Likelihood test failed: " + std::to_string(result) + " != " + std::to_string(956.0983174));
        return false;
    } else {
        logger::log_info("Likelihood test passed");
    }

    return true;
}