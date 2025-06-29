#pragma once

#include "../../../../eigen_config.h"
#include "../../types/model_element.h"
#include <vector>

class ecorr_t : public model_element_t {
private:
    
    // Per-backend/system parameters
    std::vector<string_t> backend_names_;
    std::vector<Eigen::VectorXi> backend_masks_;
    std::vector<Eigen::MatrixXd> quantization_matrices_;
    
    // Epoch information
    std::vector<std::vector<int>> epochs_per_backend_;
    int total_epochs_;
    
    // Configuration
    string_t flag_;
    double epoch_window_; // Time window for grouping TOAs into epochs (days)
    int min_toas_per_epoch_; // Minimum TOAs per epoch
    
    void detect_epochs();
    Eigen::MatrixXd create_quantization_matrix(const std::vector<int>& toa_indices) const;

public:
    
    ecorr_t();
    
    void set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json) override;
    bool is_fully_specified() const override;
    bool is_valid_parameter(const string_t& param_name) const override;
    string_t get_name() const override;
    
    void write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const override;
    
    void print_formatted(int element_index) const override;
    
    // Apply ECORR to the design matrix and coefficient priors
    void apply(const std::vector<double>& parameter_values, 
               Eigen::VectorXd& powercoeff, 
               int& start_pos, 
               double maxtspan, 
               double& uniform_prior, 
               double& freq_det) const;
    
    // Get the design matrix contribution for ECORR
    void get_design_matrix(Eigen::MatrixXd& design_matrix, int start_col) const;
    
    // Get number of ECORR coefficients
    int get_num_coefficients() const { return total_epochs_; }
};