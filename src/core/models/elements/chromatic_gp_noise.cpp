#include "chromatic_gp_noise.h"
#include "../../utils/settings.h"
#include <iostream>
#include <stdexcept>
#include <cmath>

chromatic_gp_noise_t::chromatic_gp_noise_t(const std::optional<json_node_t>& element_json)
{
    if (element_json.has_value()) {
        json_node_t json = element_json.value();
        // Check for days_per_coeff at root level (like Power Law Red Noise)
        days_per_coeff = json.get_optional_value<double>("days_per_coeff").value_or(30.0);
        
        // Also check in config section for backward compatibility
        auto config_opt = json.get_optional_value<json_node_t>("config");
        if (config_opt.has_value()) {
            days_per_coeff = config_opt.value().get_optional_value<double>("days_per_coeff").value_or(days_per_coeff);
            fixed_chromatic_idx = config_opt.value().get_optional_value<double>("fixed_chromatic_idx").value_or(4.0);
        } else {
            fixed_chromatic_idx = 4.0;  // Default to scattering-like
        }
    } else {
        days_per_coeff = 30.0;
        fixed_chromatic_idx = 4.0;  // Default to scattering-like
    }

    // num_freqs and frequencies will be calculated later in calculate_frequencies()
    num_freqs = 0;
    frequencies = Eigen::VectorXd::Zero(0);
}

void chromatic_gp_noise_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (!is_valid_parameter(name)) {
        die("Invalid parameter name for Chromatic GP Noise: " + name);
    }

    // Special handling for idx parameter with mode field
    if (name == "idx") {
        // Check mode field to determine fixed vs variable
        if (param_json.has_member("mode")) {
            std::string mode = param_json.get_value<std::string>("mode");
            
            if (mode == "fixed") {
                // Fixed idx: use the value specified in JSON, don't add to parameter_map_
                if (param_json.has_member("value")) {
                    fixed_chromatic_idx = param_json.get_value<double>("value");
                    std::cout << "  • Chromatic index: " << fixed_chromatic_idx << " (fixed)" << std::endl;
                } else {
                    throw std::invalid_argument("Fixed idx parameter (mode='fixed') must specify 'value'");
                }
                return; // Don't add to parameter_map_ since it's fixed
            } else if (mode == "variable") {
                // Variable idx: add to parameter_map_ for sampling
                std::cout << "  • Chromatic index: variable (dynamic scaling enabled)" << std::endl;
            } else {
                throw std::invalid_argument("Invalid mode for idx parameter: '" + mode + "'. Use 'fixed' or 'variable'");
            }
        } else {
            throw std::invalid_argument("idx parameter must specify 'mode' field ('fixed' or 'variable')");
        }
    }

    parameter_map_[name] = param;
}

bool chromatic_gp_noise_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {
        "log10_A",           // Log10 amplitude (Enterprise: log10_A)
        "gamma",             // Spectral index (Enterprise: gamma) 
        "idx",               // Chromatic frequency index (Enterprise: idx, default 4)
        "days_per_coeff"     // Configuration parameter
    };
    return valid_params.find(param_name) != valid_params.end();
}

bool chromatic_gp_noise_t::is_fully_specified() const
{
    // Required parameters: log10_A, gamma
    // idx is optional (can use fixed value from config)
    return get_optional_parameter("log10_A").has_value() && 
           get_optional_parameter("gamma").has_value();
}

string_t chromatic_gp_noise_t::get_name() const
{
    return "Chromatic GP Noise";
}

void chromatic_gp_noise_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const
{
    if (!is_fully_specified()) {
        return;
    }

    auto log10_A_param = get_optional_parameter("log10_A");
    auto gamma_param = get_optional_parameter("gamma");
    auto idx_param = get_optional_parameter("idx");

    if (log10_A_param && gamma_param) {
        int log10_A_index = (*log10_A_param)->get_index();
        int gamma_index = (*gamma_param)->get_index();
        
        fprintf(par_file, "# Chromatic GP Noise:\n");
        fprintf(par_file, "#   log10_A: %.6f +/- %.6f\n", parameters[log10_A_index], uncertainties[log10_A_index]);
        fprintf(par_file, "#   gamma: %.6f +/- %.6f\n", parameters[gamma_index], uncertainties[gamma_index]);
        
        if (idx_param) {
            int idx_index = (*idx_param)->get_index();
            fprintf(par_file, "#   idx: %.6f +/- %.6f\n", parameters[idx_index], uncertainties[idx_index]);
        } else {
            fprintf(par_file, "#   idx: %.6f (fixed)\n", fixed_chromatic_idx);
        }
    }
}

void chromatic_gp_noise_t::calculate_frequencies(double maxtspan)
{
    // Calculate num_freqs based on legacy approach: floor(maxtspan/days_per_coeff)
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

void chromatic_gp_noise_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const
{
    if (!is_fully_specified()) {
        return;
    }

    // Get parameters (exactly like DM noise pattern)
    const parameter_t* log10_A_param = get_parameter("log10_A");
    const parameter_t* gamma_param = get_parameter("gamma");

    double amplitude = log10_A_param->get_exp_value(parameter_values);
    double gamma = gamma_param->get_value(parameter_values);
    
    // Get chromatic index (variable or fixed)
    double chromatic_idx;
    if (is_idx_fitted()) {
        const parameter_t* idx_param = get_parameter("idx");
        chromatic_idx = idx_param->get_value(parameter_values);
    } else {
        chromatic_idx = fixed_chromatic_idx;
    }
    
    // Power law coefficients for temporal frequencies (same as DM noise)
    // NOTE: No chromatic scaling here - that goes in the design matrix per-observation
    Eigen::VectorXd chrom_coeffs = (frequencies * 365.25 / maxtspan).array().pow(-gamma);

    // Add uniform prior contribution (same as DM noise)
    if (log10_A_param->prior_type == prior_type_t::uniform) {
        uniform_prior += log(amplitude);
    }

    // Normalization (same as DM noise)
    double f1yr = 1.0 / 3.16e7;
    double pl_amp = (amplitude * amplitude / 12.0 / (M_PI * M_PI)) * std::pow(f1yr, -3) / (maxtspan * 24 * 60 * 60);

    chrom_coeffs *= pl_amp;

    // Store coefficients (same pattern as DM noise)
    powercoeff.segment(start_pos, num_freqs) += chrom_coeffs;
    powercoeff.segment(start_pos + num_freqs, num_freqs) += chrom_coeffs;

    // Frequency determinant (same pattern as DM noise)
    freq_det += 2 * powercoeff.segment(start_pos, num_freqs).array().log().sum();
    
    start_pos += 2 * num_freqs;
}

// Reconstruct chromatic GP design matrix columns when chromatic index changes
void chromatic_gp_noise_t::reconstruct_design_matrix_columns(Eigen::MatrixXd& total_matrix, 
                                                            const std::vector<double>& parameter_values,
                                                            int start_col, double maxtspan) const
{
    // Only reconstruct if chromatic index is being fitted
    if (!is_idx_fitted()) {
        return; // Fixed index, no reconstruction needed
    }

    // Get current chromatic index value
    const parameter_t* idx_param = get_parameter("idx");
    double current_chromatic_idx = idx_param->get_value(parameter_values);
    
    // Reference values for scaling calculation
    double ref_chromatic_idx = 4.0; // Reference index used in initial construction
    double ref_freq = 1400.0e6; // 1400 MHz in Hz
    
    // Calculate the scaling ratio: (new_idx / ref_idx) scaling
    // For each observation: (ref_freq/obs_freq)^current_idx / (ref_freq/obs_freq)^ref_idx 
    //                     = (ref_freq/obs_freq)^(current_idx - ref_idx)
    double idx_difference = current_chromatic_idx - ref_chromatic_idx;
    
    // Update each chromatic GP column with new scaling
    for (int i = 0; i < num_freqs; i++) {
        for (int k = 0; k < globals::pulsar->nobs; k++) {
            double obs_freq = (double)globals::pulsar->obsn[k].freqSSB;
            double time = (double)globals::pulsar->obsn[k].bat;
            
            // Calculate new chromatic scaling for this observation
            double new_chrom_scaling = std::pow(ref_freq / obs_freq, current_chromatic_idx);
            
            // Get frequency for this mode
            double freq = frequencies[i] / maxtspan;
            
            // Rebuild matrix elements with new chromatic scaling
            total_matrix(k, start_col + i) = cos(2 * M_PI * freq * time) * new_chrom_scaling;
            total_matrix(k, start_col + i + num_freqs) = sin(2 * M_PI * freq * time) * new_chrom_scaling;
        }
    }
}


