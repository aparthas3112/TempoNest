#include "temponest_v1.h"
#include "../logger.h"

double temponest_v1_t::operator()(const model_space_t& model_space, const std::vector<double>& parameter_values) const
{
    logger::log_debug("Entering TempoNest likelihood");

    double uniform_prior = 0;

    // update the residuals if we are fitting any timing model parameters
    // there is always some kind of timing model so use get_element
    auto& timing_model = model_space.get_element<timing_model_t>("Timing Model");

    timing_model.update_residuals(parameter_values);

    Eigen::VectorXd Resvec = Eigen::VectorXd::Zero(globals::pulsar->nobs);

    for (int o = 0; o < globals::pulsar->nobs; o++) {
        Resvec[o] = static_cast<double>(globals::pulsar->obsn[o].residual);
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

    // Apply EFAC if present
    if (auto efac = model_space.get_optional_element<efac_t>("EFAC")) {
        efac->get().apply(parameter_values, noise, uniform_prior);
    }

    // Apply EQUAD if present
    if (auto equad = model_space.get_optional_element<equad_t>("EQUAD")) {
        equad->get().apply(parameter_values, noise, uniform_prior);
    }

    noise = noise.array().inverse();

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
        pl_red->get().apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    // Apply DM noise if present
    if (auto pl_dm = model_space.get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        pl_dm->get().apply(parameter_values, powercoeff, start_pos, maxtspan, uniform_prior, freq_det);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get Time domain likelihood//////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double timelike = (Resvec.array().square() * noise.array()).sum();
    double tdet = -noise.array().log().sum();

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////Do Algebra/////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    logtchk("Starting algebra");

    Eigen::MatrixXd NT = TotalMatrix.array().colwise() * noise.array();

    Eigen::MatrixXd TNT = TotalMatrix.transpose() * NT;

    Eigen::VectorXd NTd = NT.transpose() * Resvec;

    logtchk("Finishing main algebra");

    if (totCoeff > 0)
        TNT.diagonal().tail(totCoeff) += powercoeff.cwiseInverse();

    // Perform Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(TNT);

    // Solve the linear system
    Eigen::VectorXd chol_solution = llt.solve(NTd);

    // Calculate log determinant
    double jointdet = 2 * llt.matrixLLT().diagonal().array().log().sum();

    double freqlike = NTd.dot(chol_solution);

    double lnewChol = -0.5 * (tdet + jointdet + freq_det + timelike - freqlike) + uniform_prior;

    logtchk("Exiting TempoNest Likelihood");

    return lnewChol;
}
