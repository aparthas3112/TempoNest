#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif

//  Copyright (C) 2013 Lindley Lentati

/*
 *    This file is part of TempoNest
 *
 *    TempoNest is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    TempoNest  is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with TempoNest.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *    If you use TempoNest and as a byproduct both Tempo2 and MultiNest
 *    then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for
 * TempoNest, Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, pp. 655-672 (bibtex:
 * 2006MNRAS.369..655H) or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4, pp. 1549-1574
 * (bibtex: 2006MNRAS.372.1549E) when discussing the timing model and MultiNest Papers here.
 */

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
#include "T2toolkit.h"
#include "TempoNest.h"
#include "eigen_config.h"
#include "namespaces/settings.h"
#include "tempo2.h"
#include "types/model.h"
#include "types/model_element.h"
using namespace std;

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context)
{

    double* DerivedParams = new double[npars];

    double result = likelihood(Cube, ndim, DerivedParams, npars, context);

    delete[] DerivedParams;

    lnew = result;
}

double likelihood(double Cube[], int ndim, double phi[], int nDerived, void* context)
{

    logtchk("Entering TempoNest likelihood");

    double uniform_prior = 0;

    int TimetoMargin = model::model_space.get_design_size();

    // update the residuals if we are fitting any timing model parameters
    // there is always some kind of timing model so use get_element
    auto& timing_model = model::get_element<timing_model_t>("Timing Model");

    timing_model.update_residuals(Cube);
    int p_count = timing_model.get_fitted_dims();

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
    if (auto efac = model::get_optional_element<efac_t>("EFAC")) {
        efac->get().apply(Cube, noise, uniform_prior, p_count);
    }

    // Apply EQUAD if present
    if (auto equad = model::get_optional_element<equad_t>("EQUAD")) {
        equad->get().apply(Cube, noise, uniform_prior, p_count);
    }

    noise = noise.array().inverse();

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Initialise TotalMatrix////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    int totalsize = model::model_space.get_total_size();

    const Eigen::MatrixXd& TotalMatrix = model::model_space.get_total_matrix();

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double maxtspan = model::model_space.get_max_tspan();
    int totCoeff = model::model_space.get_noise_size();

    Eigen::VectorXd powercoeff = Eigen::VectorXd::Zero(totCoeff);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double freq_det = 0;
    int start_pos = 0;

    // Apply red noise if present
    if (auto pl_red = model::get_optional_element<pl_red_noise_t>("Power Law Red Noise")) {
        pl_red->get().apply(Cube, powercoeff, p_count, start_pos, maxtspan, uniform_prior, freq_det);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    // Apply DM noise if present
    if (auto pl_dm = model::get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        pl_dm->get().apply(Cube, powercoeff, p_count, start_pos, maxtspan, uniform_prior, freq_det);
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
