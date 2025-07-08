#include "stochastic_solar_wind.h"
#include "../../utils/settings.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

stochastic_solar_wind_t::stochastic_solar_wind_t(const std::optional<json_node_t>& element_json)
{
    if (element_json.has_value()) {
        json_node_t json = element_json.value();
        // Check for days_per_coeff at root level
        days_per_coeff = json.get_optional_value<double>("days_per_coeff").value_or(30.0);
        
        // Also check in config section for backward compatibility
        auto config_opt = json.get_optional_value<json_node_t>("config");
        if (config_opt.has_value()) {
            days_per_coeff = config_opt.value().get_optional_value<double>("days_per_coeff").value_or(days_per_coeff);
        }
    } else {
        days_per_coeff = 30.0;
    }

    // num_freqs and frequencies will be calculated later in calculate_frequencies()
    num_freqs = 0;
    frequencies = Eigen::VectorXd::Zero(0);
}

void stochastic_solar_wind_t::calculate_frequencies(double maxtspan)
{
    // Calculate num_freqs based on time span (same as other GP models)
    num_freqs = static_cast<int>(std::floor(maxtspan / days_per_coeff));
    
    // Ensure minimum of 1 frequency
    if (num_freqs < 1) {
        num_freqs = 1;
    }
    
    // Generate frequency array: 1, 2, 3, ..., num_freqs (in units of 1/Tspan)
    frequencies = Eigen::VectorXd::Zero(num_freqs);
    for (int i = 0; i < num_freqs; i++) {
        frequencies[i] = i + 1;
    }
}

void stochastic_solar_wind_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (!is_valid_parameter(name)) {
        die("Invalid parameter name: " + name + " for " + get_name());
    }

    parameter_map_[name] = param;
}

bool stochastic_solar_wind_t::is_fully_specified() const
{
    // The element is fully specified if it has no parameters (not used) 
    // or if it has either white mode parameters or GP mode parameters
    if (parameter_map_.empty()) {
        return true;  // Element not being used
    }
    
    bool has_white = parameter_map_.find("log_amplitude") != parameter_map_.end();
    bool has_gp = (parameter_map_.find("log10_A_sw") != parameter_map_.end() && 
                   parameter_map_.find("gamma_sw") != parameter_map_.end());
    
    return has_white || has_gp;
}

bool stochastic_solar_wind_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {
        "log_amplitude",    // White noise mode
        "log10_A_sw",       // GP mode amplitude
        "gamma_sw",         // GP mode spectral index
        "days_per_coeff"    // Configuration parameter
    };
    return valid_params.find(param_name) != valid_params.end();
}

string_t stochastic_solar_wind_t::get_name() const
{
    return "Stochastic Solar Wind";
}

void stochastic_solar_wind_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const
{
    if (parameter_map_.empty()) {
        return;
    }

    // Write white mode parameters
    auto log_amplitude_param = get_optional_parameter("log_amplitude");
    if (log_amplitude_param) {
        int index = (*log_amplitude_param)->get_index();
        fprintf(par_file, "# Stochastic Solar Wind (White) - Log Amplitude: %.6f +/- %.6f\n", 
                parameters[index], uncertainties[index]);
    }
    
    // Write GP mode parameters
    auto log10_A_sw_param = get_optional_parameter("log10_A_sw");
    auto gamma_sw_param = get_optional_parameter("gamma_sw");
    
    if (log10_A_sw_param && gamma_sw_param) {
        int log10_A_index = (*log10_A_sw_param)->get_index();
        int gamma_index = (*gamma_sw_param)->get_index();
        
        fprintf(par_file, "# Stochastic Solar Wind (GP):\n");
        fprintf(par_file, "#   log10_A_sw: %.6f +/- %.6f\n", parameters[log10_A_index], uncertainties[log10_A_index]);
        fprintf(par_file, "#   gamma_sw: %.6f +/- %.6f\n", parameters[gamma_index], uncertainties[gamma_index]);
    }
}

void stochastic_solar_wind_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise) const
{
    // Main apply method for white noise mode only (called from likelihood)
    // GP mode is called separately with different signature
    if (has_white_mode()) {
        apply_white(parameter_values, noise);
    }
}

void stochastic_solar_wind_t::apply_white(const std::vector<double>& parameter_values, Eigen::VectorXd& noise) const
{
    if (parameter_map_.empty()) {
        return;  // Element not being used
    }

    auto log_amplitude_param = get_optional_parameter("log_amplitude");
    if (!log_amplitude_param) {
        return;  // Parameter not found
    }

    int param_index = (*log_amplitude_param)->get_index();
    double log_amplitude = parameter_values[param_index];
    
    // Convert from log10 space: WhiteSolarWind = 10^(log_amplitude)
    double white_solar_wind = std::pow(10.0, log_amplitude);
    
    // Get reference solar wind electron density from pulsar data
    double ne_reference = globals::pulsar->ne_sw;
    
    // Apply stochastic solar wind noise to noise vector
    // Formula: SWTerm = WhiteSolarWind × tdis2[i] / ne_sw
    // noise[i] += SWTerm² (added to inverse noise variance)
    for (int i = 0; i < globals::pulsar->nobs; i++) {
        if (globals::pulsar->obsn[i].deleted == 0) {  // Only apply to non-deleted observations
            double tdis2 = globals::pulsar->obsn[i].tdis2;
            double sw_term = white_solar_wind * tdis2 / ne_reference;
            noise[i] += sw_term * sw_term;  // Add variance component
        }
    }
}

void stochastic_solar_wind_t::apply_gp(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const
{
    if (!has_gp_mode()) {
        return;
    }

    const parameter_t* log10_A_sw_param = get_parameter("log10_A_sw");
    const parameter_t* gamma_sw_param = get_parameter("gamma_sw");

    // Get parameter values
    double log10_A_sw = log10_A_sw_param->get_value(parameter_values);
    double gamma_sw = gamma_sw_param->get_value(parameter_values);

    // Calculate temporal power spectral density (same pattern as other GP models)
    // Following TempoNest/Enterprise convention: P(f_temporal) = A² × f_temporal^(-γ)
    double amplitude = std::pow(10.0, log10_A_sw);
    
    // Power law coefficients for each temporal frequency
    Eigen::VectorXd sw_coeffs = (frequencies * 365.25 / maxtspan).array().pow(-gamma_sw);

    // Add uniform prior contribution if needed
    if (log10_A_sw_param->prior_type == prior_type_t::uniform) {
        uniform_prior += log10_A_sw; // log10_A_sw is already in log space
    }
    if (gamma_sw_param->prior_type == prior_type_t::uniform) {
        uniform_prior += gamma_sw; // gamma_sw is in linear space
    }

    // Normalization following consistent GP convention (with 12π² factor)
    double f1yr = 1.0 / 3.16e7; // 1/(1 year in seconds)
    double pl_amp = (amplitude * amplitude / 12.0 / (M_PI * M_PI)) * std::pow(f1yr, -3) / (maxtspan * 24 * 60 * 60);

    sw_coeffs *= pl_amp;

    // Store coefficients for both sine and cosine components
    // Note: The solar wind geometry and frequency scaling will be applied in the design matrix
    // during total matrix construction, similar to how DM scaling is applied
    powercoeff.segment(start_pos, num_freqs) += sw_coeffs;
    powercoeff.segment(start_pos + num_freqs, num_freqs) += sw_coeffs;

    // Update frequency determinant
    freq_det += 2 * powercoeff.segment(start_pos, num_freqs).array().log().sum();
    
    // Advance position for next model element
    start_pos += 2 * num_freqs;
}