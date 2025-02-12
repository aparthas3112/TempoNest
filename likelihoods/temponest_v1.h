
#pragma once
#include "likelihood.h"

/**
 * @brief temponest v1 likelihood implementation
 *
 * Computes log-likelihood assuming timing model and simple gaussian noise models
 */
class temponest_v1_t : public likelihood_t {
public:

    /**
     * @brief Construct temponest likelihood for given data
     */
    explicit temponest_v1_t() {}

    /**
     * @brief Compute temponest v1 log-likelihood
     *
     * Applies noise scaling from model elements and computes:
     * log(L) = -0.5 * (residual^T * Σ^-1 * residual + log|Σ|)
     * where Σ is the noise covariance matrix
     */
    double operator()(const model_space_t& model_space, const std::vector<double>& parameter_values) const override;

private:
};