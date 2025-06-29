#pragma once

#include <unordered_set>
#include <vector>
#include <optional>
#include "../../utils/settings.h"
#include "../../types/model_element.h"
#include "../../types/timing_parameter.h"

// Forward declaration for Tempo2 structures
extern "C" {
    #include "tempo2.h"
}

/**
 * Safe parameter scanner for robust timing model parameter detection.
 * Replaces unsafe fitinfo-based scanning with comprehensive validation.
 */
class parameter_scanner_t {
public:
    struct validated_parameter_t {
        param_label type;           // Parameter type (e.g., param_F0, param_JUMP)
        int index;                  // Parameter array index
        std::string label;          // Human-readable parameter name
        double value;               // Fitted parameter value
        double error;               // Fitted parameter uncertainty
        bool is_jump = false;       // True for JUMP-type parameters
        
        validated_parameter_t(param_label p, int k, const std::string& lbl, 
                            double val, double err, bool jump = false)
            : type(p), index(k), label(lbl), value(val), error(err), is_jump(jump) {}
    };

private:
    /**
     * Check if a parameter label/index combination is valid and safe to access.
     * Performs comprehensive bounds checking before array access.
     */
    bool is_valid_label_index(param_label p, int k) const;
    
    /**
     * Create a validated parameter entry with full safety checks.
     * Returns std::nullopt if parameter is invalid or should be excluded.
     */
    std::optional<validated_parameter_t> create_validated_parameter(param_label p, int k) const;
    
    /**
     * Add JUMP parameters using the separate JUMP handling system.
     * Uses psr->nJumps and psr->fitJump[] for safe iteration.
     */
    void add_jump_parameters(std::vector<validated_parameter_t>& result) const;

public:
    /**
     * Scan all fitted parameters using safe, manual parameter detection.
     * Returns validated parameters that passed all safety checks.
     * 
     * This replaces the unsafe fitinfo-based approach with a robust
     * method that works with complex .par files containing many JUMP-type parameters.
     */
    std::vector<validated_parameter_t> scan_safe_parameters() const;
};

class timing_model_t : public model_element_t {
public:

    std::vector<timing_parameter_t> parameters;
    std::vector<string_t> t2_fitted_labels;
    std::vector<bool> marginalised;
    std::vector<long double> t2_fit_values;
    std::vector<long double> t2_fit_errors;
    std::vector<std::pair<int, int>> t2_fit_indices;

    int t2_total_fit;
    int design_size;

    timing_model_t();

    // initialize all the vectors once the pulsar is loaded
    void initialise();

    // function to handle updating the residuals given the current set of parameters
    void update_residuals(const std::vector<double>& parameter_values) const;

    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;

    bool is_fully_specified() const override { return true; }

    bool is_valid_parameter(const string_t&) const override { return true; }

    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;

    string_t get_name() const override { return "Timing Model"; }

    void print_formatted(int element_index) const override
    {
        // Count how many parameters are marginalised vs included
        int marginalised_count = 0;
        int included_count = 0;
        int jump_count = 0;
        
        for (size_t i = 0; i < marginalised.size(); i++) {
            if (marginalised[i]) {
                marginalised_count++;
                // Count all JUMP-type parameters (JUMP, FDJUMP, FDJUMPDM)
                if (t2_fitted_labels[i].find("JUMP") != std::string::npos) {
                    jump_count++;
                }
            } else {
                included_count++;
            }
        }
        
        std::string extra_info;
        if (marginalised_count > 0 && included_count == 0) {
            extra_info = "Fully marginalised (" + std::to_string(marginalised_count) + " parameters)";
        } else if (marginalised_count > 0 && included_count > 0) {
            extra_info = "Partially marginalised (" + std::to_string(marginalised_count) + " marginalised, " + std::to_string(included_count) + " included)";
        } else {
            extra_info = std::to_string(included_count) + " parameters";
        }
        
        output_formatter::print_element_header(element_index, get_name(), extra_info);
        
        if (parameter_map_.empty() && marginalised_count > 0) {
            output_formatter::print_marginalised_note();
            
            // Add JUMP design philosophy if JUMP-type parameters are detected
            if (jump_count > 0) {
                std::cout << "    " << "└─ " << "\033[36m" << "Design Philosophy: JUMP-type parameters (" << jump_count 
                          << ") are always marginalized, never sampled directly" << "\033[0m" << std::endl;
            }
            return;
        }
        
        // Print included parameters in formatted way
        size_t param_count = 0;
        size_t total_params = parameter_map_.size();
        
        for (const auto& [name, param] : parameter_map_) {
            param_count++;
            bool is_last = (param_count == total_params);
            
            std::string prior_str = (param.prior_type == prior_type_t::uniform) ? "uniform" : "log-uniform";
            output_formatter::print_parameter(name, prior_str, param.min_value, param.max_value, param.include, is_last);
        }
        
        // If we have marginalised parameters, add a note
        if (marginalised_count > 0 && included_count > 0) {
            std::cout << "    " << "├─ " << marginalised_count << " additional parameters marginalised from .par file" << std::endl;
            
            // Add JUMP design philosophy if JUMP-type parameters are detected
            if (jump_count > 0) {
                std::cout << "    " << "└─ " << "\033[36m" << "Design Philosophy: JUMP-type parameters (" << jump_count 
                          << ") are always marginalized, never sampled directly" << "\033[0m" << std::endl;
            }
        }
    }
};