#include "multinest.h"
#include <json_loader.h>
#include <cfloat>
#include <filesystem>

bool multinest_settings_t::validate() const
{
    if (!sampler_settings_t::validate()) {
        return false;
    }

    return num_live > 0;
}

void multinest_settings_t::load_from_json(const json_node_t& json)
{
    // Load base class settings first
    sampler_settings_t::load_from_json(json);

    importance_sampling = json.get_optional_value<int>("importance_sampling").value_or(0);
    modal = json.get_optional_value<int>("modal").value_or(0);
    constant_efficiency = json.get_optional_value<int>("constant_efficiency").value_or(0);
    num_live = json.get_optional_value<int>("live_points").value_or(500);
    efficiency = json.get_optional_value<double>("efficiency").value_or(0.1);
    update_interval = json.get_optional_value<int>("update_interval").value_or(2000);
    num_cluster_parameters = json.get_optional_value<int>("num_cluster_parameters").value_or(1);
}

/**
 * @brief Run MultiNest nested sampling algorithm
 *
 * Sets up PolyChord configuration and runs nested sampling:
 * 1. Create output directories
 * 2. Configure global instance for callbacks
 * 3. Set up MultiNest settings (nlive, repeats, etc)
 * 4. Run sampling via C interface
 * 5. Clean up global state
 *
 * @param model Model model to sample
 */
void multinest_sampler_t::run(std::shared_ptr<model_t> model)
{
    const multinest_settings_t& mn_settings = static_cast<const multinest_settings_t&>(get_settings());

    // Create output directories
    std::filesystem::path base_dir(mn_settings.output_dir);

    try {
        std::filesystem::create_directories(base_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        die("Failed to create output directories: " + std::string(e.what()));
    }

    if (!model) {
        die("No model provided to MultiNest sampler");
    }
    model_ = model;

    sampling_parameters_ = model_->get_sampling_parameters();

    // Set the MultiNest sampling parameters

    double tol = 0.5;  // tol, defines the stopping criteria
    int ndims = model_->get_fitted_dims();

    double Ztol = -1E90;  // all the modes with logZ < Ztol are ignored
    int maxModes = 100;   // expected max no. of modes (used only for memory allocation)
    int pWrap[ndims];     // which parameters to have periodic boundary conditions?
    for (int i = 0; i < ndims; i++)
        pWrap[i] = 0;

    int seed = -1;              // random no. generator seed, if < 0 then take the seed from system clock
    int fb = 1;                 // need feedback on standard output?
    int resume = 1;             // resume from a previous job?
    int outfile = 1;            // write output files?
    int initMPI = 0;            // initialize MPI routines?, relevant only if compiling with MPI set it to F
                                // if you want your main program to handle MPI initialization
    double logZero = -DBL_MAX;  // points with loglike < logZero will be ignored by MultiNest
    int maxiter = 0;            // max no. of iterations, a non-positive value means infinity. MultiNest will
                                // terminate if either it has done max no. of iterations or convergence
                                // criterion (defined through tol) has been satisfied
    void* context = 0;          // not required by MultiNest, any additional information user wants to pass

    if (mn_settings.output_dir.size() >= 100)
        die("Root Name is too long, needs to be less than 100 characters, currently: " + std::to_string(mn_settings.output_dir.size()));

    char root[100];
    for (int r = 0; r <= mn_settings.output_dir.size(); r++) {
        root[r] = mn_settings.output_dir[r];
    }

    char* chartroot = new char[mn_settings.output_dir.length() + 1];
    std::strcpy(chartroot, mn_settings.output_dir.c_str());

    nested::run(mn_settings.importance_sampling, mn_settings.modal, mn_settings.constant_efficiency, mn_settings.num_live, tol, mn_settings.efficiency, ndims, ndims,
                mn_settings.num_cluster_parameters, maxModes, mn_settings.update_interval, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglike_wrapper, dumper, this);

    model_.reset();
}