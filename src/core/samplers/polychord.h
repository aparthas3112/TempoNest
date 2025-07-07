#pragma once

#include "sampler.h"
#include "../../external/PolyChordLite/src/polychord/interfaces.hpp"

/**
 * @brief Settings specific to the PolyChord sampling algorithm
 *
 * Contains configuration parameters specific to the PolyChord nested sampling algorithm
 */
class polychord_settings_t : public sampler_settings_t {
public:

    polychord_settings_t() = default;

    /**
     * @brief Create PolyChord settings from JSON configuration
     *
     * @param json JSON configuration node
     * @return std::unique_ptr<polychord_settings_t> Unique pointer to PolyChord settings
     */
    static std::unique_ptr<polychord_settings_t> from_json(const json_node_t& json)
    {
        auto settings = std::make_unique<polychord_settings_t>();
        settings->load_from_json(json);
        return settings;
    }

    // PolyChord specific settings

    /** @brief Number of live points (default: 1250 for high-dimensional problems) */
    int num_live = 1250;

    /** @brief Number of slice sampling steps (default: 250) */
    int num_repeats = 250;

    /** @brief Precision criterion for evidence calculation (default: 0.001) */
    double precision_criterion = 0.001;

    /** @brief Enable clustering for multimodal problems (default: true) */
    bool do_clustering = true;

    /** @brief Feedback level (0=none, 1=basic, 2=detailed) (default: 1) */
    int feedback = 1;

    /** @brief Maximum number of iterations (default: -1 for unlimited) */
    int max_ndead = -1;

    /** @brief Boost posterior weighting (default: 0.0) */
    double boost_posterior = 0.0;

    /** @brief Write posterior samples (default: true) */
    bool posteriors = true;

    /** @brief Write equally weighted samples (default: true) */
    bool equals = true;

    /** @brief Write cluster posteriors (default: true) */
    bool cluster_posteriors = true;

    /** @brief Write resume files (default: true) */
    bool write_resume = true;

    /** @brief Write parameter names file (default: true) */
    bool write_paramnames = true;

    /** @brief Read resume files if available (default: true) */
    bool read_resume = true;

    /** @brief Write statistics file (default: true) */
    bool write_stats = true;

    /** @brief Write live points file (default: false) */
    bool write_live = false;

    /** @brief Write dead points file (default: false) */
    bool write_dead = false;

    /** @brief Write prior samples (default: false) */
    bool write_prior = false;

    /** @brief Random seed (default: -1 for automatic) */
    int seed = -1;
    
    /** @brief Synchronous mode for MPI (default: false for async parallelism) */
    bool synchronous = false;

    /** @brief Enable verbose debug output for likelihood evaluations (default: false) */
    bool verbose = false;
    
    /** @brief Verbosity level as string (for JSON parsing) */
    std::string verbosity_str = "false";

    /**
     * @brief Validate PolyChord specific settings
     *
     * @return bool True if settings are valid
     */
    bool validate() const override;

    /**
     * @brief Get the sampler type name
     *
     * @return string_t The type name "polychord"
     */
    string_t get_type() const override { return "polychord"; }

protected:

    /**
     * @brief Load PolyChord specific settings from JSON
     *
     * @param json JSON configuration node
     */
    void load_from_json(const json_node_t& json) override;
};

/**
 * @brief Implementation of the PolyChord nested sampling algorithm
 *
 * PolyChord is a nested sampling algorithm particularly effective for high-dimensional
 * problems (>20 parameters). It uses slice sampling within constrained priors and
 * provides excellent scaling for complex likelihood surfaces.
 */
class polychord_sampler_t : public sampler_t {
private:

    /** @brief model being sampled */
    std::shared_ptr<model_t> model_;

    /** @brief Parameters being sampled (excluding fixed parameters) */
    std::vector<const parameter_t*> sampling_parameters_;

    /**
     * @brief Static wrapper for likelihood function to match PolyChord C interface
     *
     * @param theta Physical parameter values
     * @param nDims Number of dimensions
     * @param phi Derived parameters (unused)
     * @param nDerived Number of derived parameters (unused)
     * @return Log likelihood value
     */
    static double loglike_wrapper(double* theta, int nDims, double* phi, int nDerived);

    /**
     * @brief Static wrapper for prior function to match PolyChord C interface
     *
     * @param cube Hypercube parameter values [0,1]
     * @param theta Physical parameter values (output)
     * @param nDims Number of dimensions
     */
    static void prior_wrapper(double* cube, double* theta, int nDims);

    /**
     * @brief Static wrapper for dumper function to match PolyChord C interface
     *
     * @param nDead Number of dead points
     * @param nLive Number of live points
     * @param nPar Number of parameters
     * @param live Live points array
     * @param dead Dead points array
     * @param logweights Posterior weights
     * @param logZ Log evidence value
     * @param logZerr Error on log evidence
     */
    static void dumper_wrapper(int nDead, int nLive, int nPar, 
                              double* live, double* dead, double* logweights,
                              double logZ, double logZerr);

    /** @brief Static pointer to current sampler instance for wrapper functions */
    static polychord_sampler_t* current_instance_;

public:

    explicit polychord_sampler_t(std::unique_ptr<sampler_settings_t> settings);
    void run(std::shared_ptr<model_t> model) override;
    void output_results() override;
};