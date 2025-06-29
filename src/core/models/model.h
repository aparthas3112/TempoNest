
#pragma once
#include "../likelihood/likelihood.h"
#include "model_space.h"

/**
 * @brief Wrapper class combining model space and likelihood
 *
 * Provides unified interface for sampling and optimization
 * by combining model parameters with likelihood computation.
 */
class model_t {
public:

    /**
     * @brief Construct model likelihood
     * @param model Shared pointer to model
     * @param likelihood Shared pointer to likelihood function
     */
    model_t(std::shared_ptr<model_space_t> model, std::shared_ptr<likelihood_t> likelihood) : model_space_(std::move(model)), likelihood_(std::move(likelihood)) {}

    /**
     * @brief Get number of dimensions being sampled
     * @return Total fitted dimensions from model
     */
    int get_fitted_dims() const { return model_space_->get_fitted_dims(); }

    /**
     * @brief Compute log-likelihood for parameters
     * @param parameter_values Vector of parameter values
     * @return Log-likelihood value
     */
    double calc_loglike(const std::vector<double>& parameter_values) { return likelihood_->operator()(*model_space_, parameter_values); }

    /**
     * @brief Get parameters being sampled
     * @return Vector of pointers to non-fixed parameters
     */
    std::vector<const parameter_t*> get_sampling_parameters() const { return model_space_->get_sampling_parameters(); }

    /**
     * @brief Get read-only access to the model space
     * @return Const reference to model space
     */
    const model_space_t& get_model_space() const { return *model_space_; }

    /**
     * @brief Get read-only access to the likelihood
     * @return Const reference to likelihood
     */
    const likelihood_t& get_likelihood() const { return *likelihood_; }

private:

    std::shared_ptr<model_space_t> model_space_;
    std::shared_ptr<likelihood_t> likelihood_;
};