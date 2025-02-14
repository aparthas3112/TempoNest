#include "timing_model.h"
#include "pulsar_utils.h"

timing_model_t::timing_model_t()
{

    if (globals::pulsar == nullptr) {
        return;
    }

    initialise();
}

void timing_model_t::initialise()
{
    FitInfo* fitinfo = &(globals::pulsar->fitinfo);

    long double error_scaling = std::sqrt(globals::pulsar->fitChisq / globals::pulsar->fitNfree);

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {

        param_label p = fitinfo->paramIndex[iparam];
        const int k = fitinfo->paramCounters[iparam];

        t2_fitted_labels.push_back(globals::pulsar->param[p].shortlabel[k]);

        // be default we do marginalise over the parameters
        marginalised.push_back(true);

        // keep track of their fitted values and errors
        // if this is the phase we need to handle it explicitly
        if (p == param_ZERO) {
            t2_fit_values.push_back(globals::pulsar->offset);
            t2_fit_errors.push_back(globals::pulsar->offset_e / error_scaling);
        } else {

            long double mean = globals::pulsar->param[p].prefit[k];
            long double err = globals::pulsar->param[p].err[k] / error_scaling;

            t2_fit_values.push_back(mean);
            t2_fit_errors.push_back(err);
        }

        t2_fit_indices.push_back(std::make_pair(p, k));
    }

    // this already includes the phase offset as a parameter
    design_size = t2_fitted_labels.size();
    t2_total_fit = t2_fitted_labels.size();
}

void timing_model_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t&)
{

    timing_parameter_t timing_parameter = timing_parameter_t(param);

    auto it = std::find(t2_fitted_labels.begin(), t2_fitted_labels.end(), name);
    if (it != t2_fitted_labels.end()) {
        // Element found
        auto index = std::distance(t2_fitted_labels.begin(), it);

        // update the marginalised vector for this parameter
        marginalised[index] = false;

        // set the indices
        timing_parameter.t2_p_index = t2_fit_indices[index].first;
        timing_parameter.t2_k_index = t2_fit_indices[index].second;

        // set the long double priors for this parameter
        timing_parameter.ld_pmin = t2_fit_values[index] - 50 * t2_fit_errors[index];
        timing_parameter.ld_pmax = t2_fit_values[index] + 50 * t2_fit_errors[index];

    } else {
        // Element not found
        throw std::runtime_error("Invalid parameter name for Timing Model: " + name);
    }

    parameters.push_back(timing_parameter);
    parameter_map_[name] = param;  // Add to base class map for standard parameter access

    // every added parameters reduces the size of the design matrix by one
    design_size--;
}

void timing_model_t::update_residuals(const std::vector<double>& parameter_values) const
{

    if (parameters.size() == 0)
        return;

    for (size_t p = 0; p < parameters.size(); p++) {
        long double scalar = parameter_values[p];
        long double new_value = parameters[p].ld_pmin + (parameters[p].ld_pmax - parameters[p].ld_pmin) * scalar;

        globals::pulsar->param[parameters[p].t2_p_index].val[parameters[p].t2_k_index] = new_value;
    }

    fastformBatsAll(globals::pulsar, 1);  /* Form Barycentric arrival times */
    formResiduals(globals::pulsar, 1, 1); /* Form residuals */
}

void timing_model_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const {}
