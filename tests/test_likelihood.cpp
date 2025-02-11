#include <TempoNest.h>
#include "../logger.h"
#include "../types/model.h"
#include "tests.h"

bool run_likelihood_tests()
{
    int ndims = model::model_space.get_model_dims();
    int npars = ndims;

    double* Cube = new double[4];
    Cube[0] = 0.91;
    Cube[1] = 0.85;
    Cube[2] = 0.85;
    Cube[3] = 0.9;
    Cube[4] = 0.8;

    double* DerivedParams = new double[ndims];

    double result = likelihood(Cube, ndims, DerivedParams, npars, 0);

    delete[] DerivedParams;

    if (fabs(result - 956.0983174) > 1e-6) {
        logger::log_error("Likelihood test failed");
        return false;
    } else {
        logger::log_info("Likelihood test passed");
    }

    return true;
}