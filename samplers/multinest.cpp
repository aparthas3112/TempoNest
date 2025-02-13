#include "multinest.h"
#include <TempoNest.h>
#include <logger.h>
#include <cfloat>
#include <filesystem>
#include <fstream>
#include "../json/json_loader.h"

multinest_sampler_t::multinest_sampler_t(std::unique_ptr<sampler_settings_t> settings) : sampler_t(std::move(settings)) {}

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

void multinest_sampler_t::readtxtoutput(int ndim, std::vector<parameter_stats_t>& stats)
{
    stats.clear();
    stats = std::vector<parameter_stats_t>(ndim);
    double weightsum = 0;

    // Get filename
    std::string txt_filename = get_settings().output_dir + ".txt";

    // First pass - get means, MAP and max likelihood
    std::ifstream txt_file(txt_filename);
    if (!txt_file.is_open()) {
        throw std::runtime_error("Could not open file: " + txt_filename);
    }

    double maxlike = -1e10;
    double max_posterior = 0;

    std::string line;
    while (getline(txt_file, line)) {
        std::istringstream stream(line);
        std::vector<double> values{std::istream_iterator<double>(stream), std::istream_iterator<double>()};

        if (values.size() < ndim + 2) {
            throw std::runtime_error("Invalid line format in file");
        }

        double weight = values[0];
        double likelihood = values[1];
        weightsum += weight;

        // Update maximum likelihood point
        if (likelihood > maxlike) {
            maxlike = likelihood;
            for (int i = 0; i < ndim; i++) {
                stats[i].maximum_likelihood = values[i + 2];
            }
        }

        // Update MAP point
        if (weight > max_posterior) {
            max_posterior = weight;
            for (int i = 0; i < ndim; i++) {
                stats[i].MAP = values[i + 2];
            }
        }

        // Accumulate weighted means
        for (int i = 0; i < ndim; i++) {
            stats[i].mean += values[i + 2] * weight;
        }
    }

    // Normalize means
    for (auto& result : stats) {
        result.mean /= weightsum;
    }

    // Second pass - calculate standard deviations
    txt_file.clear();
    txt_file.seekg(0);

    while (getline(txt_file, line)) {
        std::istringstream stream(line);
        std::vector<double> values{std::istream_iterator<double>(stream), std::istream_iterator<double>()};

        double weight = values[0];
        for (int i = 0; i < ndim; i++) {
            double diff = values[i + 2] - stats[i].mean;
            stats[i].stdev += weight * diff * diff;
        }
    }

    // Finalize standard deviations
    for (auto& result : stats) {
        result.stdev = std::sqrt(result.stdev / weightsum);
    }
}

// Updates maximum likelihood if better one found in live points
void multinest_sampler_t::readphyslive(int ndim, std::vector<parameter_stats_t>& stats)
{
    std::string phys_live_filename = get_settings().output_dir + "phys_live.points";

    std::ifstream phys_live_file(phys_live_filename);
    if (!phys_live_file.is_open()) {
        throw std::runtime_error("Could not open file: " + phys_live_filename);
    }

    double maxlike = -1e10;
    std::string line;

    while (getline(phys_live_file, line)) {
        std::istringstream stream(line);
        std::vector<double> values{std::istream_iterator<double>(stream), std::istream_iterator<double>()};

        if (values.size() < ndim + 1) {  // +1 for likelihood value
            throw std::runtime_error("Invalid line format in physlive file");
        }

        double like = values[ndim];

        if (like > maxlike) {
            maxlike = like;
            for (int i = 0; i < ndim; i++) {
                stats[i].maximum_likelihood = values[i];
            }
        }
    }
}

void multinest_sampler_t::output_results()
{

    std::vector<parameter_stats_t> stats;

    int n_dims = model_->get_fitted_dims();

    readtxtoutput(n_dims, stats);
    readphyslive(n_dims, stats);

    formBatsAll(globals::pulsar, 1);       // Form Barycentric arrival times
    formResiduals(globals::pulsar, 1, 1);  // Form residuals

    TNtextOutput(globals::pulsar, 1, n_dims, get_settings().output_dir, model_, stats);

    logger::log_info("finished output");
}