#include "timing_model.h"
#include "../../utils/pulsar_utils.h"
#include "../../utils/logger.h"

timing_model_t::timing_model_t()
{

    if (globals::pulsar == nullptr) {
        return;
    }

    initialise();
}

void timing_model_t::initialise()
{
    long double error_scaling = std::sqrt(globals::pulsar->fitChisq / globals::pulsar->fitNfree);

    // Use safe parameter scanner instead of unsafe fitinfo approach
    parameter_scanner_t scanner;
    std::vector<parameter_scanner_t::validated_parameter_t> safe_parameters = scanner.scan_safe_parameters();

    // Parameter detection completed without verbose logging

    // Process all validated parameters
    for (const auto& param : safe_parameters) {
        // Store the parameter label
        t2_fitted_labels.push_back(param.label);
        
        // All parameters are marginalized by default
        marginalised.push_back(true);
        
        // Store fitted values and scaled errors
        t2_fit_values.push_back(param.value);
        t2_fit_errors.push_back(param.error / error_scaling);
        
        // Store parameter indices for later access
        t2_fit_indices.push_back(std::make_pair(param.type, param.index));
    }

    // Set design matrix size
    design_size = t2_fitted_labels.size();
    t2_total_fit = t2_fitted_labels.size();
    
    // Timing model initialization completed
}

void timing_model_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t&)
{

    // for now just dont allow fitting for jumps directly.. why would you?
    if (name.find("JUMP") != std::string::npos) {
        throw std::runtime_error("Cannot fit for jumps directly");
    }

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
        timing_parameter.ld_pmin = t2_fit_values[index] + static_cast<long double>(param.min_value) * t2_fit_errors[index];
        timing_parameter.ld_pmax = t2_fit_values[index] + static_cast<long double>(param.max_value) * t2_fit_errors[index];

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
        int index = parameters[p].get_index();
        long double scalar = parameter_values[index];
        long double new_value = parameters[p].ld_pmin + (parameters[p].ld_pmax - parameters[p].ld_pmin) * scalar;

        globals::pulsar->param[parameters[p].t2_p_index].val[parameters[p].t2_k_index] = new_value;
    }

    fastformBatsAll(globals::pulsar, 1);  /* Form Barycentric arrival times */
    formResiduals(globals::pulsar, 1, 1); /* Form residuals */
}

void timing_model_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const {}

// ============================================================================
// Parameter Scanner Implementation (Option 3: Modern C++ Safe Approach)
// ============================================================================

bool parameter_scanner_t::is_valid_label_index(param_label p, int k) const 
{
    // Sanity check: parameter type should be valid
    if (p < 0 || p >= MAX_PARAMS) {
        return false;
    }
    
    // Check if this parameter index is within bounds
    if (k < 0 || k >= globals::pulsar->param[p].aSize) {
        return false;
    }
    
    // Check if parameter is actually fitted
    if (globals::pulsar->param[p].fitFlag[k] != 1) {
        return false;
    }
    
    // Check if shortlabel exists and is valid
    if (globals::pulsar->param[p].shortlabel[k] == NULL) {
        return false;
    }
    
    return true;
}

std::optional<parameter_scanner_t::validated_parameter_t> 
parameter_scanner_t::create_validated_parameter(param_label p, int k) const 
{
    // Perform comprehensive safety checks
    if (!is_valid_label_index(p, k)) {
        return std::nullopt;
    }
    
    // Exclude DM model parameters (following legacy behavior)
    if (p == param_dmmodel) {
        return std::nullopt;
    }
    
    // Extract parameter information safely
    std::string label = std::string(globals::pulsar->param[p].shortlabel[k]);
    double value = globals::pulsar->param[p].prefit[k];
    double error = globals::pulsar->param[p].err[k];
    
    // Validate that we have reasonable values
    if (error <= 0.0) {
        logger::log_error("Parameter " + label + " has invalid error: " + std::to_string(error));
        return std::nullopt;
    }
    
    return validated_parameter_t(p, k, label, value, error, false);
}

void parameter_scanner_t::add_jump_parameters(std::vector<validated_parameter_t>& result) const 
{
    // Handle JUMP parameters separately using Tempo2's jump system
    for (int i = 0; i < globals::pulsar->nJumps; i++) {
        // Check if this jump is fitted
        if (globals::pulsar->fitJump[i] == 1) {
            std::string jump_label = "JUMP_" + std::to_string(i);
            double jump_value = globals::pulsar->jumpVal[i];
            double jump_error = globals::pulsar->jumpValErr[i];
            
            // Validate jump parameter
            if (jump_error > 0.0) {
                result.emplace_back(param_JUMP, i, jump_label, jump_value, jump_error, true);
            } else {
                logger::log_error("JUMP parameter " + jump_label + " has invalid error: " + std::to_string(jump_error));
            }
        }
    }
}

std::vector<parameter_scanner_t::validated_parameter_t> 
parameter_scanner_t::scan_safe_parameters() const 
{
    std::vector<validated_parameter_t> result;
    
    // Scan parameters without verbose logging
    
    // Manual scan through all possible parameters (like legacy TempoNest)
    for (int p = 0; p < MAX_PARAMS; p++) {
        param_label param_type = static_cast<param_label>(p);
        
        // Skip JUMP parameters - handle them separately
        if (param_type == param_JUMP) {
            continue;
        }
        
        // Scan all indices for this parameter type
        for (int k = 0; k < globals::pulsar->param[p].aSize; k++) {
            auto validated_param = create_validated_parameter(param_type, k);
            if (validated_param.has_value()) {
                result.push_back(validated_param.value());
            }
        }
    }
    
    // Handle phase offset parameter explicitly (like legacy)
    if (globals::pulsar->offset_e > 0.0) {
        result.emplace_back(param_ZERO, 0, "PHASE", globals::pulsar->offset, globals::pulsar->offset_e, false);
    }
    
    // Add JUMP parameters using separate system
    add_jump_parameters(result);
    
    // Parameter scan completed without verbose logging
    
    return result;
}
