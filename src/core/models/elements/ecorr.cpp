#include "ecorr.h"
#include "../../utils/settings.h"
#include "../../utils/output_formatter.h"
#include <iostream>
#include <algorithm>
#include <set>

ecorr_t::ecorr_t() : flag_("-sys"), epoch_window_(10.0/86400.0), total_epochs_(0), min_toas_per_epoch_(1)
{
}

void ecorr_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (name == "per_backend") {
        // Get configuration from JSON
        auto flag_opt = param_json.get_optional_value<string_t>("flag");
        flag_ = flag_opt.has_value() ? flag_opt.value() : "-sys";
        
        auto epoch_window_opt = param_json.get_optional_value<double>("epoch_window");
        epoch_window_ = epoch_window_opt.has_value() ? epoch_window_opt.value() : 10.0/86400.0; // Default: 10 seconds (legacy compatible)
        
        auto min_toas_opt = param_json.get_optional_value<int>("min_toas_per_epoch");
        min_toas_per_epoch_ = min_toas_opt.has_value() ? min_toas_opt.value() : 1; // Default: 1 TOA (legacy compatible)
        
        // Detect epochs for each backend
        detect_epochs();
        
        // Create a parameter for each backend
        for (size_t i = 0; i < backend_names_.size(); i++) {
            string_t param_name = "ecorr::" + flag_ + "::" + backend_names_[i];
            parameter_t new_param = param;
            new_param.set_index(param.get_index() + i);
            parameter_map_[param_name] = new_param;
        }
        
    } else {
        throw std::runtime_error("Invalid parameter name for ECORR: " + name);
    }
}

void ecorr_t::detect_epochs()
{
    backend_names_.clear();
    backend_masks_.clear();
    quantization_matrices_.clear();
    epochs_per_backend_.clear();
    total_epochs_ = 0;
    
    // Find all unique backend values
    std::set<string_t> unique_backends;
    std::vector<string_t> toa_backends(globals::pulsar->nobs);
    int observations_with_flag = 0;
    
    for (int o = 0; o < globals::pulsar->nobs; o++) {
        bool found = false;
        for (int f = 0; f < globals::pulsar->obsn[o].nFlags; f++) {
            string_t obs_flag(globals::pulsar->obsn[o].flagID[f]);
            if (obs_flag == flag_) {
                string_t backend_name(globals::pulsar->obsn[o].flagVal[f]);
                toa_backends[o] = backend_name;
                unique_backends.insert(backend_name);
                found = true;
                observations_with_flag++;
                break;
            }
        }
        
        if (!found) {
            toa_backends[o] = "";  // Mark as excluded
        }
    }
    
    // Validate coverage and inform user
    int observations_without_flag = globals::pulsar->nobs - observations_with_flag;
    if (observations_without_flag > 0) {
        std::cout << "ECORR " << flag_ << ": " << observations_with_flag << " TOAs included, " 
                  << observations_without_flag << " excluded (no flag)" << std::endl;
    }
    
    if (observations_with_flag == 0) {
        throw std::runtime_error("ECORR " + flag_ + ": No observations have this flag");
    }
    
    // Convert to vector and sort for consistency
    backend_names_.assign(unique_backends.begin(), unique_backends.end());
    std::sort(backend_names_.begin(), backend_names_.end());
    
    
    // For each backend, detect epochs and create quantization matrix
    for (const auto& backend : backend_names_) {
        std::vector<int> backend_toa_indices;
        
        // Find all TOAs for this backend
        for (int o = 0; o < globals::pulsar->nobs; o++) {
            if (toa_backends[o] == backend) {
                backend_toa_indices.push_back(o);
            }
        }
        
        // Create mask vector
        Eigen::VectorXi mask = Eigen::VectorXi::Zero(globals::pulsar->nobs);
        for (int idx : backend_toa_indices) {
            mask(idx) = 1;
        }
        backend_masks_.push_back(mask);
        
        // Create quantization matrix for this backend
        Eigen::MatrixXd quant_matrix = create_quantization_matrix(backend_toa_indices);
        quantization_matrices_.push_back(quant_matrix);
        
        // Store epoch information
        int num_epochs = quant_matrix.cols();
        epochs_per_backend_.push_back(std::vector<int>(num_epochs));
        total_epochs_ += num_epochs;
        
    }
    
    // Remove debug output - this information is shown in print_formatted
}

Eigen::MatrixXd ecorr_t::create_quantization_matrix(const std::vector<int>& toa_indices) const
{
    if (toa_indices.empty()) {
        return Eigen::MatrixXd::Zero(globals::pulsar->nobs, 0);
    }
    
    // Get TOA times for this backend
    std::vector<std::pair<double, int>> toa_times;
    for (int idx : toa_indices) {
        toa_times.push_back({globals::pulsar->obsn[idx].bat, idx});
    }
    
    // Sort by time
    std::sort(toa_times.begin(), toa_times.end());
    
    // Group TOAs into epochs based on time proximity
    std::vector<std::vector<int>> epochs;
    std::vector<int> current_epoch;
    current_epoch.push_back(toa_times[0].second);
    double epoch_start_time = toa_times[0].first;
    
    for (size_t i = 1; i < toa_times.size(); i++) {
        double time_diff = toa_times[i].first - epoch_start_time;
        
        if (time_diff < epoch_window_) {
            // Same epoch
            current_epoch.push_back(toa_times[i].second);
        } else {
            // New epoch
            if (current_epoch.size() >= static_cast<size_t>(min_toas_per_epoch_)) {
                epochs.push_back(current_epoch);
            }
            current_epoch.clear();
            current_epoch.push_back(toa_times[i].second);
            epoch_start_time = toa_times[i].first;
        }
    }
    
    // Add the last epoch if it has enough TOAs
    if (current_epoch.size() >= static_cast<size_t>(min_toas_per_epoch_)) {
        epochs.push_back(current_epoch);
    }
    
    // Create quantization matrix
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(globals::pulsar->nobs, epochs.size());
    
    for (size_t epoch_idx = 0; epoch_idx < epochs.size(); epoch_idx++) {
        for (int toa_idx : epochs[epoch_idx]) {
            U(toa_idx, epoch_idx) = 1.0;
        }
    }
    
    return U;
}

bool ecorr_t::is_fully_specified() const
{
    return !parameter_map_.empty();
}

bool ecorr_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"per_backend"};
    return valid_params.find(param_name) != valid_params.end();
}

string_t ecorr_t::get_name() const
{
    return "ECORR";
}

void ecorr_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& /*uncertainties*/) const
{
    for (size_t i = 0; i < backend_names_.size(); i++) {
        string_t param_name = "ecorr::" + flag_ + "::" + backend_names_[i];
        if (auto param = get_optional_parameter(param_name)) {
            // Write ECORR parameter to par file (log10 of variance)
            double log10_ecorr = parameters[param.value()->get_index()];
            fprintf(par_file, "TNECORR %s %s %g\n", flag_.c_str(), backend_names_[i].c_str(), log10_ecorr);
        }
    }
}

void ecorr_t::print_formatted(int element_index) const
{
    // Include epoch and backend info in the header
    std::string extra_info = std::to_string(backend_names_.size()) + " backends, " + 
                            std::to_string(total_epochs_) + " epochs";
    output_formatter::print_element_header(element_index, get_name(), extra_info);
    
    // Call base class to print parameters
    model_element_t::print_formatted(element_index);
}

void ecorr_t::apply(const std::vector<double>& parameter_values, 
                    Eigen::VectorXd& powercoeff, 
                    int& start_pos, 
                    double /*maxtspan*/, 
                    double& uniform_prior, 
                    double& freq_det) const
{
    if (total_epochs_ == 0) return;
    
    // Expand powercoeff if needed
    int current_size = powercoeff.size();
    int needed_size = start_pos + total_epochs_;
    if (current_size < needed_size) {
        powercoeff.conservativeResize(needed_size);
        powercoeff.tail(needed_size - current_size).setZero();
    }
    
    // Set ECORR coefficients for each epoch
    int epoch_offset = 0;
    
    for (size_t backend_idx = 0; backend_idx < backend_names_.size(); backend_idx++) {
        string_t param_name = "ecorr::" + flag_ + "::" + backend_names_[backend_idx];
        auto param = get_parameter(param_name);
        
        // Get the log10_ecorr parameter value
        double log10_ecorr = parameter_values[param->get_index()];
        double ecorr_variance = std::pow(10.0, 2.0 * log10_ecorr);
        
        // Set coefficient for each epoch of this backend
        int num_epochs_this_backend = quantization_matrices_[backend_idx].cols();
        for (int epoch_idx = 0; epoch_idx < num_epochs_this_backend; epoch_idx++) {
            powercoeff(start_pos + epoch_offset + epoch_idx) = ecorr_variance;
            freq_det += std::log(ecorr_variance);
        }
        
        // Add uniform prior term if needed
        if (param->prior_type == prior_type_t::uniform) {
            uniform_prior += std::log(std::pow(10.0, log10_ecorr));
        }
        
        epoch_offset += num_epochs_this_backend;
    }
    
    start_pos += total_epochs_;
}

void ecorr_t::get_design_matrix(Eigen::MatrixXd& design_matrix, int start_col) const
{
    if (total_epochs_ == 0) return;
    
    int current_col = start_col;
    
    // Add quantization matrices for each backend
    for (const auto& quant_matrix : quantization_matrices_) {
        int num_cols = quant_matrix.cols();
        design_matrix.block(0, current_col, design_matrix.rows(), num_cols) = quant_matrix;
        current_col += num_cols;
    }
}