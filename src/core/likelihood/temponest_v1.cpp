#include "temponest_v1.h"
#include "../utils/logger.h"
#include "../utils/settings.h"
#include "gpu_functions.h"
#include <algorithm>
#include <iostream>

double temponest_v1_t::operator()(const model_space_t& model_space, const std::vector<double>& parameter_values) const
{
    
    // NaN debugging: Track each computation step
    static int debug_call_count = 0;
    debug_call_count++;
    bool debug_this_call = globals::verbose_mode; // Only debug if verbose mode is enabled
    
    // PARAMETER LOGGING: Log parameter values for extreme cases
    if (debug_this_call) {
        logger::log_info("DEBUG STEP 0 - Parameter values [0:" + std::to_string(std::min(parameter_values.size(), size_t(10))) + "]: ");
        for (size_t i = 0; i < std::min(parameter_values.size(), size_t(10)); ++i) {
            logger::log_info("  param[" + std::to_string(i) + "] = " + std::to_string(parameter_values[i]));
        }
        if (parameter_values.size() > 10) {
            logger::log_info("  ... (total " + std::to_string(parameter_values.size()) + " parameters)");
        }
    }

    double uniform_prior = 0;

    // update the residuals if we are fitting any timing model parameters
    // there is always some kind of timing model so use get_element
    auto& timing_model = model_space.get_element<timing_model_t>("Timing Model");

    timing_model.update_residuals(parameter_values);

    Eigen::VectorXd Resvec = Eigen::VectorXd::Zero(globals::pulsar->nobs);

    for (int o = 0; o < globals::pulsar->nobs; o++) {
        Resvec[o] = static_cast<double>(globals::pulsar->obsn[o].residual);
    }
    
    if (debug_this_call) {
        double min_res = Resvec.minCoeff();
        double max_res = Resvec.maxCoeff();
        bool has_nan_res = !Resvec.allFinite();
        logger::log_info("DEBUG STEP 1 - Residuals: min=" + std::to_string(min_res) + 
                        ", max=" + std::to_string(max_res) + ", has_nan=" + (has_nan_res ? "YES" : "NO"));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get White noise vector///////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd noise = Eigen::VectorXd::Zero(globals::pulsar->nobs);

    for (int o = 0; o < globals::pulsar->nobs; o++) {
        if (!globals::use_original_errors) {
            noise[o] = globals::pulsar->obsn[o].toaErr * pow(10.0, -6);
        } else {
            noise[o] = globals::pulsar->obsn[o].origErr * pow(10.0, -6);
        }
    }
    
    if (debug_this_call) {
        double min_noise = noise.minCoeff();
        double max_noise = noise.maxCoeff();
        bool has_nan_noise = !noise.allFinite();
        logger::log_info("DEBUG STEP 2 - Initial noise: min=" + std::to_string(min_noise) + 
                        ", max=" + std::to_string(max_noise) + ", has_nan=" + (has_nan_noise ? "YES" : "NO"));
    }

    // Apply EFAC if present
    if (auto efac = model_space.get_optional_element<efac_t>("EFAC")) {
        efac->get().apply(parameter_values, noise, uniform_prior);
        if (debug_this_call) {
            double min_noise_efac = noise.minCoeff();
            double max_noise_efac = noise.maxCoeff();
            bool has_nan_efac = !noise.allFinite();
            logger::log_info("DEBUG STEP 3 - After EFAC: min=" + std::to_string(min_noise_efac) + 
                            ", max=" + std::to_string(max_noise_efac) + ", has_nan=" + (has_nan_efac ? "YES" : "NO"));
        }
    }

    // Apply Deterministic Solar Wind if present (modifies residuals)
    if (auto det_sw = model_space.get_optional_element<deterministic_solar_wind_t>("Deterministic Solar Wind")) {
        det_sw->get().apply(parameter_values, Resvec);
        if (debug_this_call) {
            double min_res_sw = Resvec.minCoeff();
            double max_res_sw = Resvec.maxCoeff();
            bool has_nan_res_sw = !Resvec.allFinite();
            logger::log_info("DEBUG STEP 4 - After Det Solar Wind: min=" + std::to_string(min_res_sw) + 
                            ", max=" + std::to_string(max_res_sw) + ", has_nan=" + (has_nan_res_sw ? "YES" : "NO"));
        }
    }

    // Apply Stochastic Solar Wind if present (adds to noise variance)
    if (auto stoch_sw = model_space.get_optional_element<stochastic_solar_wind_t>("Stochastic Solar Wind")) {
        stoch_sw->get().apply(parameter_values, noise);
        if (debug_this_call) {
            double min_noise_stoch = noise.minCoeff();
            double max_noise_stoch = noise.maxCoeff();
            bool has_nan_stoch = !noise.allFinite();
            logger::log_info("DEBUG STEP 5 - After Stoch Solar Wind: min=" + std::to_string(min_noise_stoch) + 
                            ", max=" + std::to_string(max_noise_stoch) + ", has_nan=" + (has_nan_stoch ? "YES" : "NO"));
        }
    }

    // Square the noise to get variance (noise is currently sigma, need sigma^2)
    // This matches the legacy implementation: 1/(EFAC^2 * sigma^2 + EQUAD + ...)
    noise = noise.array().square();
    
    // Apply EQUAD if present (adds to noise variance)
    if (auto equad = model_space.get_optional_element<equad_t>("EQUAD")) {
        equad->get().apply(parameter_values, noise, uniform_prior);
        if (debug_this_call) {
            double min_noise_equad = noise.minCoeff();
            double max_noise_equad = noise.maxCoeff();
            bool has_nan_equad = !noise.allFinite();
            logger::log_info("DEBUG STEP 6 - After EQUAD: min=" + std::to_string(min_noise_equad) + 
                            ", max=" + std::to_string(max_noise_equad) + ", has_nan=" + (has_nan_equad ? "YES" : "NO"));
        }
    }
    
    noise = noise.array().inverse();
    
    if (debug_this_call) {
        double min_noise_inv = noise.minCoeff();
        double max_noise_inv = noise.maxCoeff();
        bool has_nan_inv = !noise.allFinite();
        logger::log_info("DEBUG STEP 7 - After noise inverse: min=" + std::to_string(min_noise_inv) + 
                        ", max=" + std::to_string(max_noise_inv) + ", has_nan=" + (has_nan_inv ? "YES" : "NO"));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Get TotalMatrix///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    const Eigen::MatrixXd& TotalMatrix = model_space.get_total_matrix();

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double maxtspan = model_space.get_max_tspan();
    int totCoeff = model_space.get_noise_size();

    Eigen::VectorXd powercoeff = Eigen::VectorXd::Zero(totCoeff);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double freq_det = 0;
    int start_pos = 0;

    // Apply red noise if present
    if (auto pl_red = model_space.get_optional_element<pl_red_noise_t>("Power Law Red Noise")) {
        int old_start_pos = start_pos;
        double old_freq_det = freq_det;
        pl_red->get().apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
        
        if (debug_this_call) {
            int red_coeffs = start_pos - old_start_pos;
            if (red_coeffs > 0) {
                logger::log_info("DEBUG STEP 7.1 - Red Noise applied: coeffs=" + std::to_string(red_coeffs) +
                                ", freq_det_contrib=" + std::to_string(freq_det - old_freq_det));
                
                // Show some powercoeff values for red noise
                double min_red = 1e100, max_red = -1e100;
                for (int i = old_start_pos; i < start_pos; ++i) {
                    if (powercoeff[i] < min_red) min_red = powercoeff[i];
                    if (powercoeff[i] > max_red) max_red = powercoeff[i];
                }
                logger::log_info("DEBUG STEP 7.1 - Red powercoeff range: min=" + std::to_string(min_red) + 
                                ", max=" + std::to_string(max_red));
                
                // Log first few red noise coefficients
                std::string red_coeffs_str = "red_powercoeff[0:5]: ";
                for (int i = old_start_pos; i < std::min(start_pos, old_start_pos + 5); ++i) {
                    red_coeffs_str += std::to_string(powercoeff[i]) + " ";
                }
                logger::log_info("DEBUG STEP 7.1 - " + red_coeffs_str);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    // Apply DM noise if present
    if (auto pl_dm = model_space.get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        int old_start_pos = start_pos;
        double old_freq_det = freq_det;
        pl_dm->get().apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
        
        if (debug_this_call) {
            int dm_coeffs = start_pos - old_start_pos;
            if (dm_coeffs > 0) {
                logger::log_info("DEBUG STEP 7.2 - DM Noise applied: coeffs=" + std::to_string(dm_coeffs) +
                                ", freq_det_contrib=" + std::to_string(freq_det - old_freq_det));
                
                // Show some powercoeff values for DM noise
                double min_dm = 1e100, max_dm = -1e100;
                for (int i = old_start_pos; i < start_pos; ++i) {
                    if (powercoeff[i] < min_dm) min_dm = powercoeff[i];
                    if (powercoeff[i] > max_dm) max_dm = powercoeff[i];
                }
                logger::log_info("DEBUG STEP 7.2 - DM powercoeff range: min=" + std::to_string(min_dm) + 
                                ", max=" + std::to_string(max_dm));
                
                // Log first few DM noise coefficients  
                std::string dm_coeffs_str = "dm_powercoeff[0:5]: ";
                for (int i = old_start_pos; i < std::min(start_pos, old_start_pos + 5); ++i) {
                    dm_coeffs_str += std::to_string(powercoeff[i]) + " ";
                }
                logger::log_info("DEBUG STEP 7.2 - " + dm_coeffs_str);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Chromatic GP Noise//////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    // Apply Chromatic GP Noise if present  
    if (auto chrom_gp = model_space.get_optional_element<chromatic_gp_noise_t>("Chromatic GP Noise")) {
        int old_start_pos = start_pos;
        double old_freq_det = freq_det;
        chrom_gp->get().apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
        
        if (debug_this_call) {
            int chrom_coeffs = start_pos - old_start_pos;
            if (chrom_coeffs > 0) {
                logger::log_info("DEBUG STEP 7.3 - Chromatic GP Noise applied: coeffs=" + std::to_string(chrom_coeffs) +
                                ", freq_det_contrib=" + std::to_string(freq_det - old_freq_det));
                
                // Show some powercoeff values for chromatic GP noise
                double min_chrom = 1e100, max_chrom = -1e100;
                for (int i = old_start_pos; i < start_pos; ++i) {
                    if (powercoeff[i] < min_chrom) min_chrom = powercoeff[i];
                    if (powercoeff[i] > max_chrom) max_chrom = powercoeff[i];
                }
                logger::log_info("DEBUG STEP 7.3 - Chromatic GP powercoeff range: min=" + std::to_string(min_chrom) + 
                                ", max=" + std::to_string(max_chrom));
                
                // Log first few chromatic GP noise coefficients
                std::string chrom_coeffs_str = "chrom_gp_powercoeff[0:5]: ";
                for (int i = old_start_pos; i < std::min(start_pos, old_start_pos + 5); ++i) {
                    chrom_coeffs_str += std::to_string(powercoeff[i]) + " ";
                }
                logger::log_info("DEBUG STEP 7.3 - " + chrom_coeffs_str);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Stochastic Solar Wind GP//////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    // Apply Stochastic Solar Wind GP if present
    if (auto stoch_sw = model_space.get_optional_element<stochastic_solar_wind_t>("Stochastic Solar Wind")) {
        if (stoch_sw->get().has_gp_mode()) {
            int old_start_pos = start_pos;
            double old_freq_det = freq_det;
            stoch_sw->get().apply_gp(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
            
            if (debug_this_call) {
                int sw_coeffs = start_pos - old_start_pos;
                if (sw_coeffs > 0) {
                    logger::log_info("DEBUG STEP 7.4 - Solar Wind GP applied: coeffs=" + std::to_string(sw_coeffs) +
                                    ", freq_det_contrib=" + std::to_string(freq_det - old_freq_det));
                    
                    // Show some powercoeff values for solar wind GP
                    double min_sw = 1e100, max_sw = -1e100;
                    for (int i = old_start_pos; i < start_pos; ++i) {
                        if (powercoeff[i] < min_sw) min_sw = powercoeff[i];
                        if (powercoeff[i] > max_sw) max_sw = powercoeff[i];
                    }
                    logger::log_info("DEBUG STEP 7.4 - Solar Wind GP powercoeff range: min=" + std::to_string(min_sw) + 
                                    ", max=" + std::to_string(max_sw));
                    
                    // Log first few solar wind GP coefficients
                    std::string sw_coeffs_str = "sw_gp_powercoeff[0:5]: ";
                    for (int i = old_start_pos; i < std::min(start_pos, old_start_pos + 5); ++i) {
                        sw_coeffs_str += std::to_string(powercoeff[i]) + " ";
                    }
                    logger::log_info("DEBUG STEP 7.4 - " + sw_coeffs_str);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////ECORR Epochs///////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    // Apply ECORR if present - handle multiple instances
    int total_ecorr_coeffs = 0;
    double ecorr_freq_det_start = freq_det;
    const auto& elements = model_space.get_elements();
    
    for (const auto& [name, element] : elements) {
        if (name.find("ECORR") == 0) {
            auto* ecorr = element->as<ecorr_t>();
            if (ecorr) {
                int old_start_pos = start_pos;
                ecorr->apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
                total_ecorr_coeffs += (start_pos - old_start_pos);
            }
        }
    }
    
    if (debug_this_call && total_ecorr_coeffs > 0) {
        logger::log_info("DEBUG STEP 7.3 - ECORR applied: coeffs=" + std::to_string(total_ecorr_coeffs) +
                        ", freq_det_contrib=" + std::to_string(freq_det - ecorr_freq_det_start));
        
        // Show some powercoeff values for ECORR
        double min_ecorr = 1e100, max_ecorr = -1e100;
        int ecorr_start = start_pos - total_ecorr_coeffs;
        for (int i = ecorr_start; i < start_pos; ++i) {
            if (powercoeff[i] < min_ecorr) min_ecorr = powercoeff[i];
            if (powercoeff[i] > max_ecorr) max_ecorr = powercoeff[i];
        }
        logger::log_info("DEBUG STEP 7.3 - ECORR powercoeff range: min=" + std::to_string(min_ecorr) + 
                        ", max=" + std::to_string(max_ecorr));
    }
    
    // COMPREHENSIVE POWERCOEFF SUMMARY
    if (debug_this_call && totCoeff > 0) {
        double overall_min = 1e100, overall_max = -1e100;
        int extreme_count = 0;
        std::vector<int> extreme_indices;
        
        for (int i = 0; i < totCoeff; ++i) {
            if (powercoeff[i] < overall_min) overall_min = powercoeff[i];
            if (powercoeff[i] > overall_max) overall_max = powercoeff[i];
            
            // Track extremely large values (> 1e15)
            if (powercoeff[i] > 1e15) {
                extreme_count++;
                if (extreme_indices.size() < 5) extreme_indices.push_back(i);
            }
        }
        
        logger::log_info("DEBUG STEP 7.4 - POWERCOEFF SUMMARY: total_coeffs=" + std::to_string(totCoeff) +
                        ", overall_min=" + std::to_string(overall_min) + ", overall_max=" + std::to_string(overall_max));
        logger::log_info("DEBUG STEP 7.4 - Extreme values (>1e15): count=" + std::to_string(extreme_count));
        
        if (!extreme_indices.empty()) {
            std::string extreme_details = "Extreme coefficient indices: ";
            for (size_t i = 0; i < extreme_indices.size(); ++i) {
                int idx = extreme_indices[i];
                extreme_details += "idx[" + std::to_string(idx) + "]=" + std::to_string(powercoeff[idx]) + " ";
            }
            logger::log_info("DEBUG STEP 7.4 - " + extreme_details);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get Time domain likelihood//////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double timelike = (Resvec.array().square() * noise.array()).sum();
    double tdet = -noise.array().log().sum();
    
    if (debug_this_call) {
        bool timelike_finite = std::isfinite(timelike);
        bool tdet_finite = std::isfinite(tdet);
        logger::log_info("DEBUG STEP 8 - Time domain: timelike=" + std::to_string(timelike) + 
                        " (finite=" + (timelike_finite ? "YES" : "NO") + "), tdet=" + std::to_string(tdet) +
                        " (finite=" + (tdet_finite ? "YES" : "NO") + ")");
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////Do Algebra/////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // Reconstruct chromatic GP matrix for variable chromatic index
    model_space.reconstruct_chromatic_gp_matrix(parameter_values);
    
    logtchk("Starting algebra");

    double likelihood = 0.0;
    if (globals::use_gpu) {
        if (debug_this_call) {
            logger::log_info("DEBUG STEP 8.0 - Using GPU path, optimized=" + std::string(globals::use_gpu_optimized ? "YES" : "NO"));
        }
        
        if (globals::use_gpu_optimized) {
            // Initialize enhanced GPU cache if not already done
            if (!gpu_data::isInitialized()) {
                gpu_data::initialize_with_cache(TotalMatrix, TotalMatrix.rows() + totCoeff);
            }
            likelihood = performAlgebraWithArrayFireGPU_optimized(noise, Resvec, powercoeff, totCoeff, tdet, freq_det, timelike, uniform_prior);
        } else {
            likelihood = performAlgebraWithArrayFireGPU(TotalMatrix, noise, Resvec, powercoeff, totCoeff, tdet, freq_det, timelike, uniform_prior);
        }
        
        if (debug_this_call) {
            bool gpu_likelihood_finite = std::isfinite(likelihood);
            logger::log_info("DEBUG STEP 8.4 - GPU result: likelihood=" + std::to_string(likelihood) + 
                            " (finite=" + std::string(gpu_likelihood_finite ? "YES" : "NO") + ")");
        }
    } else {
        if (debug_this_call) {
            logger::log_info("DEBUG STEP 8.0 - Using CPU path");
        }
        Eigen::MatrixXd NT = TotalMatrix.array().colwise() * noise.array();

        Eigen::MatrixXd TNT = TotalMatrix.transpose() * NT;

        Eigen::VectorXd NTd = NT.transpose() * Resvec;

        logtchk("Finishing main algebra");

        if (totCoeff > 0) {
            // Check powercoeff for NaN/zero before inversion
            if (debug_this_call) {
                double min_powercoeff = powercoeff.minCoeff();
                double max_powercoeff = powercoeff.maxCoeff();
                bool has_nan_powercoeff = !powercoeff.allFinite();
                logger::log_info("DEBUG STEP 8.1 - Powercoeff: min=" + std::to_string(min_powercoeff) + 
                                ", max=" + std::to_string(max_powercoeff) + ", has_nan=" + std::string(has_nan_powercoeff ? "YES" : "NO"));
            }
            
            TNT.diagonal().tail(totCoeff) += powercoeff.cwiseInverse();
        }

        // Perform Cholesky decomposition
        Eigen::LLT<Eigen::MatrixXd> llt(TNT);
        
        if (debug_this_call) {
            bool chol_success = (llt.info() == Eigen::Success);
            logger::log_info("DEBUG STEP 8.2 - Cholesky decomposition: success=" + std::string(chol_success ? "YES" : "NO"));
        }

        // Solve the linear system
        Eigen::VectorXd chol_solution = llt.solve(NTd);

        // Calculate log determinant
        double jointdet = 2 * llt.matrixLLT().diagonal().array().log().sum();

        double freqlike = NTd.dot(chol_solution);
        
        if (debug_this_call) {
            bool jointdet_finite = std::isfinite(jointdet);
            bool freqlike_finite = std::isfinite(freqlike);
            logger::log_info("DEBUG STEP 8.3 - Matrix ops: jointdet=" + std::to_string(jointdet) + 
                            " (finite=" + std::string(jointdet_finite ? "YES" : "NO") + "), freqlike=" + std::to_string(freqlike) +
                            " (finite=" + std::string(freqlike_finite ? "YES" : "NO") + ")");
        }

        likelihood = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;
    }
    
    if (debug_this_call) {
        bool likelihood_finite = std::isfinite(likelihood);
        logger::log_info("DEBUG STEP 9 - Final likelihood: " + std::to_string(likelihood) + 
                        " (finite=" + (likelihood_finite ? "YES" : "NO") + ")");
        logger::log_info("DEBUG COMPONENTS - tdet=" + std::to_string(tdet) + 
                        ", freq_det=" + std::to_string(freq_det) + 
                        ", timelike=" + std::to_string(timelike) + 
                        ", uniform_prior=" + std::to_string(uniform_prior));
    }
    
    logtchk("Exiting TempoNest Likelihood");

    // Record the end of likelihood call for profiling
    
    return likelihood;
}
