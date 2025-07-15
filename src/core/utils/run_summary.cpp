#include "run_summary.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <ctime>
#include <regex>
#include "settings.h"
#include "logger.h"
#include "../samplers/multinest.h"
#include "../models/elements/timing_model.h"
#include "../models/elements/power_law_red_noise.h"
#include "../models/elements/power_law_dm_noise.h"
#include "../types/parameter.h"
#include "../likelihood/gpu_functions.h"

run_summary_t::run_summary_t(const std::string& output_directory) 
    : output_dir_(output_directory), timing_started_(false) {
}

void run_summary_t::start_timing() {
    start_time_ = std::chrono::steady_clock::now();
    timing_started_ = true;
}

void run_summary_t::finish_and_write_summary(std::shared_ptr<model_t> model, 
                                            const sampler_settings_t& sampler_settings) {
    if (!timing_started_) {
        logger::log_error("Timing was never started - cannot write accurate summary");
        return;
    }
    
    end_time_ = std::chrono::steady_clock::now();
    
    // Calculate runtime in seconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_);
    double runtime_seconds = duration.count() / 1000.0;
    
    // Generate summary content
    std::string summary_content = generate_summary_content(model, sampler_settings, runtime_seconds);
    
    // Extract directory from output_root (MultiNest uses it as a filename prefix)
    std::string output_dir = std::filesystem::path(sampler_settings.output_root).parent_path().string();
    if (output_dir.empty()) {
        output_dir = ".";  // Current directory if no path specified
    }
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir);
    
    // Write summary file in the same directory as chains
    std::string summary_filepath = output_dir + "/run_summary.txt";
    std::ofstream summary_file(summary_filepath);
    
    if (summary_file.is_open()) {
        summary_file << summary_content;
        summary_file.close();
        logger::log_info("Run summary written to: " + summary_filepath);
    } else {
        logger::log_error("Failed to write run summary to: " + summary_filepath);
    }
    
    // Also write parameter names file
    write_parameter_names(model, sampler_settings);
}

std::string run_summary_t::generate_summary_content(std::shared_ptr<model_t> model,
                                                   const sampler_settings_t& sampler_settings,
                                                   double runtime_seconds) const {
    std::stringstream ss;
    
    // Header
    ss << "┌──────────────────────────────────────────────────────────────────────┐\n";
    ss << "│                        TEMPONEST RUN SUMMARY                         │\n";
    ss << "└──────────────────────────────────────────────────────────────────────┘\n\n";
    
    // Basic run information
    ss << "Run Information:\n";
    ss << "  • Completion Time: " << get_timestamp() << "\n";
    ss << "  • Total Runtime: " << format_runtime(runtime_seconds) << "\n";
    ss << "  • Output Directory: " << output_dir_ << "\n\n";
    
    // Global settings
    ss << "Global Settings:\n";
    ss << "  • Debug Mode: " << (globals::debug ? "ON" : "OFF") << "\n";
    ss << "  • Tempo2 Iterations: " << globals::num_tempo2_its << "\n";
    ss << "  • Original Errors: " << (globals::use_original_errors ? "ON" : "OFF") << "\n";
    ss << "  • Test Mode: " << (globals::test_mode ? "ON" : "OFF") << "\n";
#ifdef HAVE_ARRAYFIRE
    ss << "  • GPU Acceleration: " << (globals::use_gpu ? "ON" : "OFF") << "\n";
#else
    ss << "  • GPU Acceleration: NOT AVAILABLE\n";
#endif
    ss << "\n";
    
    // Sampler configuration
    ss << "Sampler Configuration:\n";
    ss << "  • Type: " << sampler_settings.get_type() << "\n";
    ss << "  • Sampling: " << (sampler_settings.sample ? "ENABLED" : "DISABLED") << "\n";
    ss << "  • Output Root: " << sampler_settings.output_root << "\n";
    
    // Add specific sampler settings if available
    if (auto* mn_settings = dynamic_cast<const multinest_settings_t*>(&sampler_settings)) {
        ss << "  • Live Points: " << mn_settings->num_live << "\n";
        ss << "  • Efficiency: " << mn_settings->efficiency << "\n";
        ss << "  • Update Interval: " << mn_settings->update_interval << "\n";
    }
    ss << "\n";
    
    // Frequency coefficient summary
    ss << "Frequency Coefficient Summary:\n";
    ss << "  • Time Span: " << std::fixed << std::setprecision(1) << model->get_model_space().get_max_tspan() << " days\n";
    
    const auto& elements = model->get_model_space().get_elements();
    for (const auto& [name, element] : elements) {
        if (element->get_name() == "Power Law Red Noise") {
            if (auto* red_noise = element->as<pl_red_noise_t>()) {
                ss << "  • Red Noise: " << red_noise->num_freqs << " frequencies (" 
                   << red_noise->days_per_coeff << " days/coeff, "
                   << (2 * red_noise->num_freqs) << " total coeffs)\n";
            }
        } else if (element->get_name() == "Power Law DM Noise") {
            if (auto* dm_noise = element->as<pl_dm_noise_t>()) {
                ss << "  • DM Noise: " << dm_noise->num_freqs << " frequencies (" 
                   << dm_noise->days_per_coeff << " days/coeff, "
                   << (2 * dm_noise->num_freqs) << " total coeffs)\n";
            }
        }
    }
    ss << "\n";
    
    // Model information
    ss << "┌──────────────────────────────────────────────────────────────────────┐\n";
    ss << "│                            MODEL ELEMENTS                            │\n";
    ss << "└──────────────────────────────────────────────────────────────────────┘\n\n";
    
    size_t i = 0;
    for (const auto& [name, element] : elements) {
        i++;
        ss << "[" << i << "] " << element->get_name();
        
        // Add element-specific information
        if (element->get_name() == "Power Law Red Noise") {
            if (auto* red_noise = element->as<pl_red_noise_t>()) {
                ss << " (" << red_noise->num_freqs << " frequencies, " 
                   << red_noise->days_per_coeff << " days/coeff, "
                   << (2 * red_noise->num_freqs) << " total coeffs)";
            } else {
                ss << " (Red noise process)";
            }
        } else if (element->get_name() == "Power Law DM Noise") {
            if (auto* dm_noise = element->as<pl_dm_noise_t>()) {
                ss << " (" << dm_noise->num_freqs << " frequencies, " 
                   << dm_noise->days_per_coeff << " days/coeff, "
                   << (2 * dm_noise->num_freqs) << " total coeffs)";
            } else {
                ss << " (DM noise process)";
            }
        } else if (element->get_name() == "EFAC") {
            ss << " (Error scaling)";
        } else if (element->get_name() == "EQUAD") {
            ss << " (Quadrature error)";
        } else if (element->get_name() == "Timing Model") {
            // Get timing model info
            auto fitted_dims = element->get_fitted_dims();
            
            // Try to cast to timing_model_t to get marginalisation info
            if (auto* timing_model = element->as<timing_model_t>()) {
                int marginalised_count = 0;
                for (bool is_marginalised : timing_model->marginalised) {
                    if (is_marginalised) marginalised_count++;
                }
                
                int total_dims = marginalised_count + fitted_dims;
                
                if (fitted_dims == 0 && total_dims > 0) {
                    ss << " (Fully marginalised, " << total_dims << " parameters)";
                } else if (marginalised_count > 0 && fitted_dims > 0) {
                    ss << " (Partially marginalised, " << marginalised_count 
                       << " marginalised, " << fitted_dims << " included)";
                } else {
                    ss << " (" << fitted_dims << " parameters)";
                }
            } else {
                ss << " (" << fitted_dims << " parameters)";
            }
        }
        ss << "\n";
        
        // List parameters for this element
        const auto& params = element->get_parameters();
        for (size_t j = 0; j < params.size(); ++j) {
            const auto* param = params[j];
            std::string connector = (j == params.size() - 1) ? "└─" : "├─";
            
            ss << "    " << connector << " " << param->name << ": ";
            
            if (param->include) {
                std::string prior_str = (param->prior_type == prior_type_t::uniform) ? "uniform" : "log-uniform";
                ss << prior_str << " [" << param->min_value << ", " << param->max_value << "] ✓ included";
            } else {
                ss << "✗ excluded";
            }
            ss << "\n";
        }
        
        if (i < elements.size() - 1) {
            ss << "\n";
        }
    }
    
    // Summary statistics
    ss << "\n";
    ss << "Model Summary:\n";
    ss << "  • Total Elements: " << elements.size() << "\n";
    ss << "  • Fitted Parameters: " << model->get_fitted_dims() << "\n";
    ss << "  • Total Parameters: " << model->get_fitted_dims() << "\n";
    
    // Data information
    ss << "\nData Information:\n";
    ss << "  • Pulsar: " << globals::pulsar->name << "\n";
    ss << "  • Observations: " << globals::pulsar->nobs << "\n";
    ss << "  • Time Span: " << std::fixed << std::setprecision(1) 
       << (globals::pulsar->param[param_finish].val[0] - globals::pulsar->param[param_start].val[0])/365.25 
       << " years\n";
    
    ss << "\n";
    ss << "Generated by TempoNest " << get_timestamp() << "\n";
    
    return ss.str();
}

std::string run_summary_t::get_timestamp() const {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    
    std::stringstream ss;
    ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

std::string run_summary_t::format_runtime(double seconds) const {
    if (seconds < 60.0) {
        return std::to_string(static_cast<int>(seconds)) + " seconds";
    } else if (seconds < 3600.0) {
        int minutes = static_cast<int>(seconds / 60.0);
        int remaining_seconds = static_cast<int>(seconds) % 60;
        return std::to_string(minutes) + "m " + std::to_string(remaining_seconds) + "s";
    } else {
        int hours = static_cast<int>(seconds / 3600.0);
        int remaining_minutes = static_cast<int>((seconds - hours * 3600) / 60.0);
        return std::to_string(hours) + "h " + std::to_string(remaining_minutes) + "m";
    }
}

void run_summary_t::write_parameter_names(std::shared_ptr<model_t> model,
                                         const sampler_settings_t& sampler_settings) const {
    // Extract directory from output_root
    std::string output_dir = std::filesystem::path(sampler_settings.output_root).parent_path().string();
    if (output_dir.empty()) {
        output_dir = ".";  // Current directory if no path specified
    }
    
    // Get sampling parameters in order
    auto parameters = model->get_sampling_parameters();
    
    // Write parameter names file
    std::string paramnames_filepath = output_dir + "/paramnames.txt";
    std::ofstream paramnames_file(paramnames_filepath);
    
    if (paramnames_file.is_open()) {
        for (const auto* param : parameters) {
            std::string param_name;
            
            if (!param->id.empty()) {
                // Use the full parameter ID which includes context
                param_name = param->id;
                
                // Make human-readable: convert "per_flag::-fe::0" to "efac_fe_0" 
                if (param_name.find("per_flag::") == 0) {
                    std::string element_name = param->get_parent()->get_name();
                    if (element_name == "EFAC") {
                        param_name = std::regex_replace(param_name, std::regex("per_flag::([^:]+)::(.+)"), "efac_$1_$2");
                        param_name = std::regex_replace(param_name, std::regex("-"), "");  // Remove dashes
                    } else if (element_name == "EQUAD") {
                        param_name = std::regex_replace(param_name, std::regex("per_flag::([^:]+)::(.+)"), "equad_$1_$2");
                        param_name = std::regex_replace(param_name, std::regex("-"), "");  // Remove dashes
                    }
                } else {
                    // For simple parameter names, add element-specific prefixes
                    std::string element_name = param->get_parent()->get_name();
                    if (element_name == "Power Law Red Noise") {
                        param_name = "red_" + param_name;
                    } else if (element_name == "Power Law DM Noise") {
                        param_name = "dm_" + param_name;
                    } else if (element_name == "Chromatic GP Noise") {
                        if (param->name == "log10_A") {
                            param_name = "chrom_gp_amp";
                        } else if (param->name == "idx") {
                            param_name = "chrom_gp_index";
                        } else if (param->name == "gamma") {
                            param_name = "chrom_gp_gamma";
                        } else {
                            param_name = "chrom_gp_" + param_name;
                        }
                    }
                }
            } else {
                // For simple parameters, add model element prefix
                std::string element_name = param->get_parent()->get_name();
                if (element_name == "Power Law Red Noise") {
                    param_name = "red_" + param->name;
                } else if (element_name == "Power Law DM Noise") {
                    param_name = "dm_" + param->name;
                } else if (element_name == "Chromatic GP Noise") {
                    if (param->name == "log10_A") {
                        param_name = "chrom_gp_amp";
                    } else if (param->name == "idx") {
                        param_name = "chrom_gp_index";
                    } else if (param->name == "gamma") {
                        param_name = "chrom_gp_gamma";
                    } else {
                        param_name = "chrom_gp_" + param->name;
                    }
                } else if (element_name == "EFAC") {
                    param_name = "efac_" + param->name;
                } else if (element_name == "EQUAD") {
                    param_name = "equad_" + param->name;
                } else if (element_name == "Timing Model") {
                    param_name = param->name;  // Keep timing parameter names as-is
                } else {
                    param_name = param->name;  // Fallback to original name
                }
            }
            
            paramnames_file << param_name << "\n";
        }
        paramnames_file.close();
        logger::log_info("Parameter names written to: " + paramnames_filepath);
    } else {
        logger::log_error("Failed to write parameter names to: " + paramnames_filepath);
    }
}