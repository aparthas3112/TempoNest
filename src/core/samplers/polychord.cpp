#include "polychord.h"
#include "../utils/logger.h"
#include "../utils/settings.h"
#include "../models/model.h"
#include <algorithm>
#include <filesystem>
#include <thread>
#include <atomic>
#include <chrono>
#include <cmath>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

// Static member initialization
polychord_sampler_t* polychord_sampler_t::current_instance_ = nullptr;

// PolyChord settings implementation

bool polychord_settings_t::validate() const
{
    // Call base class validation first
    if (!sampler_settings_t::validate()) {
        return false;
    }

    // Validate PolyChord specific settings
    if (num_live <= 0) {
        logger::log_error("PolyChord num_live must be positive");
        return false;
    }

    if (num_repeats <= 0) {
        logger::log_error("PolyChord num_repeats must be positive");
        return false;
    }

    if (precision_criterion <= 0.0) {
        logger::log_error("PolyChord precision_criterion must be positive");
        return false;
    }

    if (feedback < 0 || feedback > 2) {
        logger::log_error("PolyChord feedback must be 0, 1, or 2");
        return false;
    }

    return true;
}

void polychord_settings_t::load_from_json(const json_node_t& json)
{
    // Load base class settings first
    sampler_settings_t::load_from_json(json);

    // Load PolyChord specific settings with optimized defaults for pulsar timing
    num_live = json.get_optional_value<int>("num_live").value_or(200);
    num_repeats = json.get_optional_value<int>("num_repeats").value_or(25);
    precision_criterion = json.get_optional_value<double>("precision_criterion").value_or(0.001);
    do_clustering = json.get_optional_value<bool>("do_clustering").value_or(false);  // Optimized for pulsar timing
    feedback = json.get_optional_value<int>("feedback").value_or(1);
    max_ndead = json.get_optional_value<int>("max_ndead").value_or(-1);
    boost_posterior = json.get_optional_value<double>("boost_posterior").value_or(5.0);  // Critical for pulsar timing
    posteriors = json.get_optional_value<bool>("posteriors").value_or(true);
    equals = json.get_optional_value<bool>("equals").value_or(true);
    cluster_posteriors = json.get_optional_value<bool>("cluster_posteriors").value_or(false);  // Match do_clustering
    write_resume = json.get_optional_value<bool>("write_resume").value_or(true);
    write_paramnames = json.get_optional_value<bool>("write_paramnames").value_or(true);
    read_resume = json.get_optional_value<bool>("read_resume").value_or(true);
    write_stats = json.get_optional_value<bool>("write_stats").value_or(true);
    write_live = json.get_optional_value<bool>("write_live").value_or(false);
    write_dead = json.get_optional_value<bool>("write_dead").value_or(false);
    write_prior = json.get_optional_value<bool>("write_prior").value_or(false);
    seed = json.get_optional_value<int>("seed").value_or(-1);
    synchronous = json.get_optional_value<bool>("synchronous").value_or(false);
    
    // Handle verbose as either bool or string
    if (auto verbose_bool = json.get_optional_value<bool>("verbose")) {
        verbose = verbose_bool.value();
        verbosity_str = verbose ? "true" : "false";
    } else if (auto verbose_string = json.get_optional_value<std::string>("verbose")) {
        verbosity_str = verbose_string.value();
        // Convert string to bool for backward compatibility
        verbose = (verbosity_str == "true" || verbosity_str == "full");
    } else {
        // Default to normal verbosity
        verbose = true;
        verbosity_str = "true";
    }

    // Alternative parameter names for backwards compatibility
    if (auto live_points = json.get_optional_value<int>("live_points")) {
        num_live = live_points.value();
    }
    if (auto nlive = json.get_optional_value<int>("nlive")) {
        num_live = nlive.value();
    }
    if (auto precision = json.get_optional_value<double>("precision")) {
        precision_criterion = precision.value();
    }
    if (auto clustering = json.get_optional_value<bool>("clustering")) {
        do_clustering = clustering.value();
    }
}

// PolyChord sampler implementation

polychord_sampler_t::polychord_sampler_t(std::unique_ptr<sampler_settings_t> settings)
    : sampler_t(std::move(settings))
{
}

void polychord_sampler_t::run(std::shared_ptr<model_t> model)
{
    model_ = model;
    sampling_parameters_ = model_->get_sampling_parameters();
    
    // Set static instance pointer for wrapper functions
    current_instance_ = this;

    const auto& settings = static_cast<const polychord_settings_t&>(get_settings());
    
    // Set global verbosity level based on sampler settings
    if (settings.verbosity_str == "false") {
        globals::set_verbosity_level(globals::VerbosityLevel::QUIET);
    } else if (settings.verbosity_str == "true") {
        globals::set_verbosity_level(globals::VerbosityLevel::NORMAL);
    } else if (settings.verbosity_str == "full") {
        globals::set_verbosity_level(globals::VerbosityLevel::FULL);
    } else {
        // Default to normal for any other value
        globals::set_verbosity_level(globals::VerbosityLevel::NORMAL);
    }

    // Create output directories
    if (!create_output_directories()) {
        die("Failed to create output directories for PolyChord");
    }

    // Create additional PolyChord-specific directories
    std::filesystem::path output_path(settings.output_root);
    std::string base_dir = output_path.parent_path().string() + "/";
    
    try {
        std::filesystem::create_directories(base_dir + "clusters");
        std::filesystem::create_directories(base_dir + "chains");
        logger::log_info("Created PolyChord subdirectories");
    } catch (const std::filesystem::filesystem_error& e) {
        logger::log_error("Failed to create PolyChord subdirectories: " + std::string(e.what()));
        die("Failed to create PolyChord subdirectories");
    }

    logger::log_info("Starting PolyChord nested sampling");
    logger::log_info("Number of parameters: " + std::to_string(sampling_parameters_.size()));

    // Configure PolyChord settings
    Settings pc_settings;
    pc_settings.nDims = static_cast<int>(sampling_parameters_.size());
    pc_settings.nDerived = 0;  // No derived parameters for now
    pc_settings.nlive = settings.num_live;
    
    // Use user-specified num_repeats from JSON configuration
    pc_settings.num_repeats = settings.num_repeats;
    
    pc_settings.precision_criterion = settings.precision_criterion;
    
    // FIX 2: Critical settings from legacy implementation
    pc_settings.do_clustering = false;  // Always false for pulsar timing (legacy)
    pc_settings.boost_posterior = 5.0;  // CRITICAL: Legacy value that worked
    
    pc_settings.feedback = settings.feedback;
    pc_settings.max_ndead = settings.max_ndead;
    pc_settings.posteriors = settings.posteriors;
    pc_settings.equals = settings.equals;
    pc_settings.cluster_posteriors = false;  // Disabled in legacy
    pc_settings.write_resume = settings.write_resume;
    pc_settings.write_paramnames = settings.write_paramnames;
    pc_settings.read_resume = settings.read_resume;
    pc_settings.write_stats = settings.write_stats;
    pc_settings.write_live = settings.write_live;
    pc_settings.write_dead = settings.write_dead;
    pc_settings.write_prior = settings.write_prior;
    pc_settings.seed = settings.seed;
    
    // Set critical missing parameters with proper defaults
    pc_settings.logzero = -1e30;
    pc_settings.nprior = -1;
    pc_settings.nfail = -1;
    
    // FIX 3: Use synchronous mode for stability (but respect user setting)
    pc_settings.synchronous = settings.synchronous;
    
    // Additional critical parameter for MPI builds
    pc_settings.maximise = false;
    
    // FIX 4: Use legacy compression factor (proven value)
    pc_settings.compression_factor = 0.36787944117144233;
    
    // Set grade parameters (required for slice sampling)
    pc_settings.grade_frac = {1.0};  // Single grade with all parameters
    pc_settings.grade_dims = {pc_settings.nDims};  // All dimensions in one grade

    // Log final settings with legacy compatibility notes
    logger::log_info("Number of live points: " + std::to_string(pc_settings.nlive));
    logger::log_info("Number of repeats (user-specified): " + std::to_string(pc_settings.num_repeats));
    logger::log_info("Precision criterion: " + std::to_string(pc_settings.precision_criterion));
    logger::log_info("Boost posterior (legacy value): " + std::to_string(pc_settings.boost_posterior));
    logger::log_info("Clustering disabled: " + std::string(pc_settings.do_clustering ? "false" : "true"));
    logger::log_info("Compression factor (legacy): " + std::to_string(pc_settings.compression_factor));

    // Set output files
    pc_settings.base_dir = base_dir;
    pc_settings.file_root = output_path.filename().string();

    // Ensure base_dir ends with separator
    if (!pc_settings.base_dir.empty() && pc_settings.base_dir.back() != '/') {
        pc_settings.base_dir += "/";
    }

    logger::log_info("Output directory: " + pc_settings.base_dir);
    logger::log_info("File root: " + pc_settings.file_root);

    // FIX 5: Proper MPI handling and enhanced error detection
    #ifdef HAVE_MPI
    // Initialize MPI if not already done
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        logger::log_info("Initializing MPI for PolyChord");
        MPI_Init(nullptr, nullptr);
    }
    
    // Create MPI communicator
    MPI_Comm polychord_comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &polychord_comm);
    logger::log_info("Using MPI communicator for PolyChord");
    #endif
    
    // Run PolyChord with comprehensive error checking
    try {
        logger::log_info("Starting PolyChord nested sampling...");
        logger::log_info("Expected output directory: " + pc_settings.base_dir);
        logger::log_info("Expected file prefix: " + pc_settings.file_root);
        
        #ifdef HAVE_MPI
        run_polychord(loglike_wrapper, prior_wrapper, dumper_wrapper, pc_settings, polychord_comm);
        MPI_Comm_free(&polychord_comm);
        #else
        run_polychord(loglike_wrapper, prior_wrapper, dumper_wrapper, pc_settings);
        #endif
        
        // FIX 6: Verify output files were created
        std::string stats_file = pc_settings.base_dir + pc_settings.file_root + ".stats";
        std::string resume_file = pc_settings.base_dir + pc_settings.file_root + ".resume";
        
        if (!std::filesystem::exists(stats_file)) {
            throw std::runtime_error("PolyChord failed to create stats file: " + stats_file + 
                                   "\nPossible causes: 1) Likelihood function errors, 2) MPI issues, 3) Directory permissions");
        }
        
        logger::log_info("PolyChord sampling completed successfully");
        logger::log_info("Output files created in: " + pc_settings.base_dir);
        
    } catch (const std::exception& e) {
        logger::log_error("PolyChord sampling failed: " + std::string(e.what()));
        logger::log_error("Debug info: nDims=" + std::to_string(pc_settings.nDims) + 
                         ", nlive=" + std::to_string(pc_settings.nlive) + 
                         ", num_repeats=" + std::to_string(pc_settings.num_repeats));
        logger::log_error("Check: 1) Likelihood function implementation, 2) Parameter bounds, 3) Output directory permissions");
        throw;
    } catch (...) {
        logger::log_error("PolyChord sampling failed with unknown error");
        logger::log_error("This may indicate: 1) Memory issues, 2) Library linking problems, 3) MPI communication errors");
        throw;
    }

    // Clear static instance pointer
    current_instance_ = nullptr;
}

void polychord_sampler_t::output_results()
{
    const auto& settings = static_cast<const polychord_settings_t&>(get_settings());
    
    logger::log_info("PolyChord results written to: " + settings.output_root);
    
    // Additional result processing could be added here
    // For now, PolyChord handles its own output file generation
}

// Static wrapper functions

double polychord_sampler_t::loglike_wrapper(double* theta, int nDims, double* /*phi*/, int /*nDerived*/)
{
    static int call_count = 0;
    call_count++;
    
    if (!current_instance_) {
        logger::log_error("PolyChord loglike_wrapper called without valid instance");
        return -1e30;
    }

    try {
        // Convert parameter array to vector
        std::vector<double> params(theta, theta + nDims);
        
        // DEBUG: Log some likelihood evaluations (only for FULL verbosity)
        static int debug_count = 0;
        if (globals::verbosity_level == globals::VerbosityLevel::FULL && 
            (debug_count < 5 || debug_count % 1000 == 0)) {
            std::string param_str = "";
            for (int i = 0; i < std::min(5, nDims); ++i) {
                param_str += std::to_string(params[i]) + " ";
            }
            logger::log_info("PolyChord likelihood eval #" + std::to_string(debug_count) + 
                            ": params[0:4]=" + param_str);
        }
        
        // Calculate log likelihood using the model
        double loglike = current_instance_->model_->calc_loglike(params);
        
        // DEBUG: Log likelihood values (only for FULL verbosity)
        if (globals::verbosity_level == globals::VerbosityLevel::FULL && 
            (debug_count < 5 || debug_count % 1000 == 0)) {
            logger::log_info("PolyChord likelihood result #" + std::to_string(debug_count) + 
                            ": loglike=" + std::to_string(loglike));
        }
        debug_count++;
        
        return loglike;
    } catch (const std::exception& e) {
        logger::log_error("Error in PolyChord loglike_wrapper: " + std::string(e.what()));
        return -1e30;
    } catch (...) {
        logger::log_error("Unknown error in PolyChord loglike_wrapper");
        return -1e30;
    }
}

void polychord_sampler_t::prior_wrapper(double* cube, double* theta, int nDims)
{
    if (!current_instance_) {
        logger::log_error("PolyChord prior_wrapper called without valid instance");
        return;
    }

    try {
        const auto& parameters = current_instance_->sampling_parameters_;
        
        // DEBUG: Log parameter transformation details on first call (only for FULL verbosity)
        static bool first_call = true;
        static int prior_call_count = 0;
        if (globals::verbosity_level == globals::VerbosityLevel::FULL && first_call) {
            logger::log_info("PolyChord prior_wrapper: nDims=" + std::to_string(nDims) + 
                            ", parameters.size()=" + std::to_string(parameters.size()));
            for (int i = 0; i < nDims && i < static_cast<int>(parameters.size()); ++i) {
                logger::log_info("Parameter " + std::to_string(i) + ": " + parameters[i]->name + 
                                " [" + std::to_string(parameters[i]->min_value) + 
                                ", " + std::to_string(parameters[i]->max_value) + "]");
            }
            first_call = false;
        }
        
        // DEBUG: Log some transformations (only for FULL verbosity)
        if (globals::verbosity_level == globals::VerbosityLevel::FULL && prior_call_count < 3) {
            std::string transform_str = "";
            for (int i = 0; i < std::min(5, nDims); ++i) {
                if (i < static_cast<int>(parameters.size())) {
                    double min_val = parameters[i]->min_value;
                    double max_val = parameters[i]->max_value;
                    double physical = min_val + cube[i] * (max_val - min_val);
                    transform_str += std::to_string(cube[i]) + "→" + std::to_string(physical) + " ";
                }
            }
            logger::log_info("PolyChord prior transform #" + std::to_string(prior_call_count) + 
                            ": cube→physical: " + transform_str);
        }
        prior_call_count++;
        
        // Transform from unit hypercube [0,1] to physical parameter space
        for (int i = 0; i < nDims; ++i) {
            if (i < static_cast<int>(parameters.size())) {
                double min_val = parameters[i]->min_value;
                double max_val = parameters[i]->max_value;
                theta[i] = min_val + cube[i] * (max_val - min_val);
            } else {
                // Safety fallback
                logger::log_error("PolyChord prior_wrapper: parameter index " + std::to_string(i) + 
                                " out of range (size=" + std::to_string(parameters.size()) + ")");
                theta[i] = cube[i];
            }
        }
    } catch (const std::exception& e) {
        logger::log_error("Error in PolyChord prior_wrapper: " + std::string(e.what()));
    } catch (...) {
        logger::log_error("Unknown error in PolyChord prior_wrapper");
    }
}

void polychord_sampler_t::dumper_wrapper(int nDead, int /*nLive*/, int /*nPar*/, 
                                       double* /*live*/, double* /*dead*/, double* /*logweights*/,
                                       double logZ, double logZerr)
{
    // This function is called periodically during PolyChord sampling
    // Can be used for runtime monitoring, progress updates, etc.
    
    if (current_instance_) {
        // Progress updates based on verbosity level
        if (globals::verbosity_level >= globals::VerbosityLevel::NORMAL) {
            // For NORMAL verbosity: show progress every 100 points
            if (nDead > 0 && nDead % 100 == 0) {
                logger::log_info("PolyChord progress: " + std::to_string(nDead) + 
                                " dead points, logZ = " + std::to_string(logZ) + 
                                " ± " + std::to_string(logZerr));
            }
            
            // Also log the first point to confirm sampling started
            if (nDead == 1) {
                logger::log_info("PolyChord sampling started - first dead point processed");
            }
        }
        
        // More frequent updates for FULL verbosity
        if (globals::verbosity_level == globals::VerbosityLevel::FULL) {
            if (nDead > 0 && nDead % 10 == 0) {
                logger::log_info("PolyChord detailed progress: " + std::to_string(nDead) + 
                                " dead points, logZ = " + std::to_string(logZ) + 
                                " ± " + std::to_string(logZerr));
            }
            
            // Log early points for debugging
            if (nDead <= 5 && nDead > 1) {
                logger::log_info("PolyChord early progress: " + std::to_string(nDead) + 
                                " dead points, logZ = " + std::to_string(logZ));
            }
        }
    }
}