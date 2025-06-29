#pragma once
#include <memory>
#include <vector>
#include "../models/model_space.h"

/**
 * @brief Abstract base class for likelihood functions
 *
 * Provides interface for computing log-likelihood values given a model
 * and parameter values.
 */
class likelihood_t {
public:

    virtual ~likelihood_t() = default;

    /**
     * @brief Compute log-likelihood
     * @param model Model containing elements and parameters
     * @param parameters Vector of parameter values in physical space
     * @return Log-likelihood value
     */
    virtual double operator()(const model_space_t& model_space, const std::vector<double>& parameters) const = 0;

protected:

    likelihood_t() = default;
};
