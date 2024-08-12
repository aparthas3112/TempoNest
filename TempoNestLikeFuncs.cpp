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

void* globalcontext;

void assigncontext(void* context)
{
    globalcontext = context;
}

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context)
{

    double* DerivedParams = new double[npars];

    double result = NewLRedMarginLogLike(Cube, ndim, DerivedParams, npars, context);

    delete[] DerivedParams;

    lnew = result;
}

double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived, void* context)
{

    logtchk("Entering TempoNest likelihood");

    double uniformpriorterm = 0;

    int TimetoMargin = model::design_size;

    Eigen::VectorXd Resvec = Eigen::VectorXd::Zero(model::pulsar->nobs);
    double EQUAD = 0;
    double EFAC = 1;

    int pcount = 0;

    for (int o = 0; o < model::pulsar->nobs; o++) {

        Resvec[o] = (double)model::pulsar->obsn[o].residual;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get White Noise vector///////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (model::efac.has_value()) {
        efac_element* efac = model::efac.value()->as<efac_element>();
        EFAC = efac->global.value().get_exp_value(Cube[pcount]);

        if (efac->global.value().prior_type == prior_type_t::uniform) {
            uniformpriorterm += log(EFAC);
        }

        pcount++;
    }

    // printf("Equad %i \n", ((MNStruct *)globalcontext)->numFitEQUAD);
    if (model::equad.has_value()) {
        equad_element* equad = model::equad.value()->as<equad_element>();
        double value = equad->global.value().get_exp_value(Cube[pcount]);
        value = value * value;

        EQUAD = value;
        if (equad->global.value().prior_type == prior_type_t::uniform) {
            uniformpriorterm += 0.5 * log(value);
        }

        pcount++;
    }

    Eigen::VectorXd Noise = Eigen::VectorXd::Zero(model::pulsar->nobs);

    for (int o = 0; o < model::pulsar->nobs; o++) {
        double EFACterm = 0;
        double noiseval = 0;

        if (!settings::use_original_errors) {
            noiseval = model::pulsar->obsn[o].toaErr;
        } else {
            noiseval = model::pulsar->obsn[o].origErr;
        }

        EFACterm = (noiseval * pow(10.0, -6)) * EFAC;

        Noise[o] = 1.0 / (pow(EFACterm, 2) + EQUAD);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Initialise TotalMatrix////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    int totalsize = model::total_size;

    Eigen::MatrixXd TotalMatrix = model::total_matrix;

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double maxtspan = model::max_tspan;
    int totCoeff = model::noise_size;

    Eigen::VectorXd powercoeff = Eigen::VectorXd::Zero(totCoeff);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd freqs = Eigen::VectorXd::Zero(totCoeff);
    Eigen::VectorXd DMVec = Eigen::VectorXd::Zero(model::pulsar->nobs);

    double freqdet = 0;
    int startpos = 0;

    if (model::pl_red_noise.has_value()) {

        pl_red_noise_element* pl = model::pl_red_noise.value()->as<pl_red_noise_element>();

        const int halfCoeff = pl->num_freqs;

        // Calculate and assign frequencies
        freqs.segment(startpos, halfCoeff) = pl->frequencies / maxtspan;

        freqs.segment(startpos + halfCoeff, halfCoeff) = freqs.segment(startpos, halfCoeff);

        double Tspan = maxtspan;
        double f1yr = 1.0 / 3.16e7;

        double redamp = pl->amplitude.get_exp_value(Cube[pcount++]);
        double redindex = pl->spectral_index.get_value(Cube[pcount++]);

        if (pl->amplitude.prior_type == prior_type_t::uniform) {
            uniformpriorterm += log(redamp);
        }

        for (int i = 0; i < halfCoeff; i++) {

            double rho = (redamp * redamp / 12.0 / (M_PI * M_PI)) * pow(f1yr, (-3)) *
                         pow(freqs[i] * 365.25, (-redindex)) / (Tspan * 24 * 60 * 60);

            powercoeff[i] += rho;
            powercoeff[i + halfCoeff] += rho;
        }

        startpos = 2 * pl->num_freqs;

        freqdet += 2 * powercoeff.segment(0, halfCoeff).array().log().sum();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (model::pl_dm_noise.has_value()) {

        pl_dm_noise_element* pl = model::pl_dm_noise.value()->as<pl_dm_noise_element>();

        const int halfCoeff = pl->num_freqs;

        double DMKappa = 2.410 * pow(10.0, -16);

        for (int o = 0; o < model::pulsar->nobs; o++) {
            DMVec[o] = 1.0 / (DMKappa * pow((double)model::pulsar->obsn[o].freqSSB, 2));
        }

        // Calculate and assign frequencies
        freqs.segment(startpos, halfCoeff) = pl->frequencies / maxtspan;

        freqs.segment(startpos + halfCoeff, halfCoeff) = freqs.segment(startpos, halfCoeff);

        double DMamp = pl->amplitude.get_exp_value(Cube[pcount++]);
        double DMindex = pl->spectral_index.get_value(Cube[pcount++]);

        double f1yr = 1.0 / 3.16e7;

        if (pl->amplitude.prior_type == prior_type_t::uniform) {
            uniformpriorterm += log(DMamp);
        }

        for (int i = 0; i < halfCoeff; i++) {

            double rho = (DMamp * DMamp) * pow(f1yr, (-3)) *
                         pow(freqs[startpos + i] * 365.25, (-DMindex)) / (maxtspan * 24 * 60 * 60);
            powercoeff[startpos + i] += rho;
            powercoeff[startpos + i + halfCoeff] += rho;
        }

        freqdet += 2 * powercoeff.segment(startpos, halfCoeff).array().log().sum();

        startpos += 2 * pl->num_freqs;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get Time domain likelihood//////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double timelike = (Resvec.array().square() * Noise.array()).sum();
    double tdet = -Noise.array().log().sum();

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////Do Algebra/////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    logtchk("Starting algebra");

    Eigen::MatrixXd NT = TotalMatrix.array().colwise() * Noise.array();

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

    double lnewChol = -0.5 * (tdet + jointdet + freqdet + timelike - freqlike) + uniformpriorterm;

    logtchk("Exiting TempoNest Likelihood");

    return lnewChol;
}
