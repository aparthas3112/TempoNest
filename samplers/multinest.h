#pragma once

#include "multinest_interface.h"
#include "sampler.h"

/**
 * @brief Settings specific to the MultiNest sampling algorithm
 *
 * Contains configuration parameters specific to the MultiNest nested sampling algorithm
 */
class multinest_settings_t : public sampler_settings_t {
public:

    multinest_settings_t() = default;

    /**
     * @brief Create MultiNest settings from JSON configuration
     *
     * @param json JSON configuration node
     * @return std::unique_ptr<multinest_settings_t> Unique pointer to MultiNest settings
     */
    static std::unique_ptr<multinest_settings_t> from_json(const json_node_t& json)
    {
        auto settings = std::make_unique<multinest_settings_t>();
        settings->load_from_json(json);
        return settings;
    }

    // IS: flag to use importance sampling, will be supported by upcoming release of multinest, at
    // which point it should be set to 1
    int importance_sampling;

    // modal: flag to allow multinest to search for multiple modes in the data. 1 = multimodal, 0 =
    // single mode
    int modal;

    // ceff: flag to set multinest to constant efficiency mode.  Adjusts sampling to maintain the
    // efficiency set by the efr parameter. This is usefull for large dimensional problems (>
    // 20dim), however the accuracy of the evidence suffers if importance sampling isn't used.
    int constant_efficiency;

    // nlive: the number of live points used by multinest.  THe more you have the more it explores
    // the parameter space, but the longer it takes to sample. As a guide depending on the
    // dimensionality of the problem (this is the numer of parameters sampled after discounting
    // those to be marginalised over): Less than 5 dimensions: 100 Between 5 and 10 dimensions: 200
    // Between 10 and 50: 500
    // More than 50: 1000
    int num_live;

    // efr: The sampling efficiency.  The Lower this number the more carefully multinest explores
    // the parameter space, so the longer it takes. The value depends on the dimensionality and the
    // goal. For parameter estimation it can be set higher than if an accurate Evidence value is
    // required. As a rough guide: Less than 10 dimensions: 0.8 (parameter estimation), 0.3
    // (Evidence evaluation) Between 10 and 20 dimensions: 0.3 (parameter estimation), 0.1 (Evidence
    // evaluation) Between 20 and 50: 0.1 (parameter estimation), 0.05 (Evidence evaluation) More
    // than 50: 0.05 (parameter estimation), 0.01 (Evidence evaluation) These numbers may well be
    // adjusted with experience, for parameter estimation of a single modal "blob", these may be
    // much lower than required. Also: In constant efficiency mode this must be set lower, to 0.05
    // for D<50 and 0.01 D>50.
    double efficiency;

    // do sampling? 0 = No, 1 = Yes
    bool sample;

    // updInt: How often the output files are updated
    int update_interval;

    // nClsPar:  Number of parameters to cluster over when doing multi modal analysis
    int num_cluster_parameters;

    /**
     * @brief Validate MultiNest specific settings
     *
     * @return bool True if settings are valid
     */
    bool validate() const override;

protected:

    /**
     * @brief Load MultiNest specific settings from JSON
     *
     * @param json JSON configuration node
     */
    void load_from_json(const json_node_t& json) override;
};

/**
 * @brief Implementation of the MultiNest nested sampling algorithm
 *
 * MultiNest is a nested sampling algorithm that uses slice sampling
 * to draw points from within iso-likelihood contours. This implementation
 * provides C++ wrappers around the MultiNest Fortran library.
 */
class multinest_sampler_t : public sampler_t {
private:

    /** @brief model being sampled */
    std::shared_ptr<model_t> model_;

    /** @brief Parameters being sampled (excluding fixed parameters) */
    std::vector<const parameter_t*> sampling_parameters_;

    /**
     * @brief Static wrapper for likelihood function to match MultiNest C interface
     *
     * @param Cube Hypercube parameter values
     * @param ndim Number of dimensions
     * @param phi Derived parameters (unused)
     * @param npars Number of derived parameters (unused)
     */
    static void loglike_wrapper(double* Cube, int& ndim, int& npars, double& lnew, void* context)
    {

        // Cast the context back to the object instance and call the instance method
        multinest_sampler_t* self = static_cast<multinest_sampler_t*>(context);
        std::vector<const parameter_t*> parameters = self->model_->get_sampling_parameters();

        std::vector<double> params(ndim);

        for (int i = 0; i < ndim; ++i) {
            double physical_value = parameters[i]->min_value + Cube[i] * (parameters[i]->max_value - parameters[i]->min_value);
            Cube[i] = physical_value;
            params[i] = physical_value;
        }

        lnew = self->model_->calc_loglike(params);
    }

    /************************************************* dumper routine
     * ******************************************************/

    // The dumper routine will be called every updInt*10 iterations
    // MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way
    // they want
    //
    //
    // Arguments:
    //
    // nSamples 						= total number of samples in posterior distribution
    // nlive 						= total number of live points
    // nPar 						= total number of parameters (free + derived)
    // physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points
    // (physical parameters plus derived parameters) along with their loglikelihood values
    // posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points.
    // Each sample has nPar parameters (physical + derived) along with the their loglike value &
    // posterior probability paramConstr[1][4*nPar]: paramConstr[0][0] to paramConstr[0][nPar - 1] 	=
    // mean values of the parameters paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard
    // deviation of the parameters paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit
    // (maxlike) parameters paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP
    // (maximum-a-posteriori) parameters maxLogLike						= maximum loglikelihood value
    // logZ							= log evidence value
    // logZerr						= error on log evidence value
    // context						void pointer, any additional information

    static void dumper(int& nSamples, int& nlive, int& nPar, double** physLive, double** posterior, double** paramConstr, double& maxLogLike, double& logZ, double& logZerr, void* context) {}

    // functions for reading output
    void readphyslive(int ndim, std::vector<parameter_stats_t>& stats);
    void readtxtoutput(int ndim, std::vector<parameter_stats_t>& stats);

public:

    explicit multinest_sampler_t(std::unique_ptr<sampler_settings_t> settings);
    void run(std::shared_ptr<model_t> model) override;
    void output_results() override;
};
