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
#include <iterator>
#include <sstream>
#include <vector>
#include "T2toolkit.h"
#include "TempoNest.h"
#include "eigen_config.h"
#include "tempo2.h"

#ifdef HAVE_MLAPACK
    #include <mpblas_mpfr.h>
    #include <mplapack_mpfr.h>
#endif
#include <iostream>

using namespace std;

void* globalcontext;

// void SmallNelderMeadOptimum(int nParameters, double *ML);
// void  WriteProfileFreqEvo(std::string longname, int &ndim, int profiledimstart);

void assigncontext(void* context)
{
    globalcontext = context;
}

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context)
{

    for (int p = 0; p < ndim; p++) {

        Cube[p] = (((MNStruct*)globalcontext)->PriorsArray[p + ndim] -
                   ((MNStruct*)globalcontext)->PriorsArray[p]) *
                      Cube[p] +
                  ((MNStruct*)globalcontext)->PriorsArray[p];
    }

    double* DerivedParams = new double[npars];

    // double result = NewLRedMarginLogLike(ndim, Cube, npars, DerivedParams, context);
    double result = NewLRedMarginLogLike(Cube, ndim, DerivedParams, npars, context);

    delete[] DerivedParams;

    lnew = result;
}

// double  NewLRedMarginLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void
// *context){
double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived, void* context)
{

    logtchk("Entering TempoNest likelihood");

    double uniformpriorterm = 0;
    clock_t startClock, endClock;

    int TimetoMargin = ((MNStruct*)globalcontext)->TimetoMargin;

    Eigen::VectorXd Resvec = Eigen::VectorXd::Zero(((MNStruct*)globalcontext)->pulse->nobs);
    Eigen::VectorXd EQUAD = Eigen::VectorXd::Zero(((MNStruct*)globalcontext)->systemcount);
    Eigen::VectorXd EFAC = Eigen::VectorXd::Ones(((MNStruct*)globalcontext)->systemcount);

    int pcount = 0;

    for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {

        Resvec[o] = (double)((MNStruct*)globalcontext)->pulse->obsn[o].residual;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get White Noise vector///////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)globalcontext)->numFitEFAC == 1) {

        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
            EFAC(o) = pow(10.0, Cube[pcount]);
            if (((MNStruct*)globalcontext)->EFACPriorType == 1) {
                uniformpriorterm += log(EFAC(o));
            }
        }
        pcount++;

    } else if (((MNStruct*)globalcontext)->numFitEFAC > 1) {

        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
            EFAC(o) = pow(10.0, Cube[pcount]);
            if (((MNStruct*)globalcontext)->EFACPriorType == 1) {
                uniformpriorterm += log(EFAC(o));
            }
            pcount++;
        }
    }

    // printf("Equad %i \n", ((MNStruct *)globalcontext)->numFitEQUAD);
    if (((MNStruct*)globalcontext)->numFitEQUAD == 1) {
        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
            EQUAD(o) = pow(10.0, 2 * Cube[pcount]);
            if (((MNStruct*)globalcontext)->EQUADPriorType == 1) {
                uniformpriorterm += log(pow(10.0, Cube[pcount]));
            }
        }
        pcount++;
    } else if (((MNStruct*)globalcontext)->numFitEQUAD > 1) {
        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {

            if (((MNStruct*)globalcontext)->includeEQsys[o] == 1) {
                // printf("Cube: %i %i %g \n", o, pcount, Cube[pcount]);
                EQUAD(o) = pow(10.0, 2 * Cube[pcount]);
                if (((MNStruct*)globalcontext)->EQUADPriorType == 1) {
                    uniformpriorterm += log(pow(10.0, Cube[pcount]));
                }
                pcount++;
            }
        }
    }

    Eigen::VectorXd Noise = Eigen::VectorXd::Zero(((MNStruct*)globalcontext)->pulse->nobs);

    double DMKappa = 2.410 * pow(10.0, -16);

    for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {
        double EFACterm = 0;
        double noiseval = 0;

        if (((MNStruct*)globalcontext)->useOriginalErrors == 0) {
            noiseval = ((MNStruct*)globalcontext)->pulse->obsn[o].toaErr;
        } else if (((MNStruct*)globalcontext)->useOriginalErrors == 1) {
            noiseval = ((MNStruct*)globalcontext)->pulse->obsn[o].origErr;
        }

        EFACterm = (noiseval * pow(10.0, -6)) * EFAC(((MNStruct*)globalcontext)->sysFlags[o]);

        Noise[o] = 1.0 / (pow(EFACterm, 2) + EQUAD(((MNStruct*)globalcontext)->sysFlags[o]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Initialise TotalMatrix////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    int totalsize = ((MNStruct*)globalcontext)->totalsize;

    Eigen::MatrixXd TotalMatrix(((MNStruct*)context)->pulse->nobs, ((MNStruct*)context)->totalsize);

    for (int i = 0; i < ((MNStruct*)globalcontext)->pulse->nobs; i++) {
        for (int j = 0; j < ((MNStruct*)globalcontext)->totalsize; j++) {
            TotalMatrix(i, j) = ((MNStruct*)globalcontext)
                                    ->StoredTMatrix[i + j * ((MNStruct*)context)->pulse->nobs];
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double maxtspan = ((MNStruct*)globalcontext)->Tspan;
    double averageTSamp = 2 * maxtspan / ((MNStruct*)globalcontext)->pulse->nobs;

    int FitRedCoeff = 2 * (((MNStruct*)globalcontext)->numFitRedCoeff);
    int FitDMCoeff = 2 * (((MNStruct*)globalcontext)->numFitDMCoeff);

    int totCoeff = ((MNStruct*)globalcontext)->totCoeff;

    Eigen::VectorXd powercoeff = Eigen::VectorXd::Zero(totCoeff);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd freqs = Eigen::VectorXd::Zero(totCoeff);
    Eigen::VectorXd DMVec = Eigen::VectorXd::Zero(((MNStruct*)globalcontext)->pulse->nobs);

    double freqdet = 0;
    int startpos = 0;

    if (((MNStruct*)globalcontext)->incRED > 0) {

        for (int i = 0; i < FitRedCoeff / 2; i++) {

            freqs[startpos + i] = (double)((MNStruct*)globalcontext)->sampleFreq[i] / maxtspan;
            freqs[startpos + i + FitRedCoeff / 2] = freqs[startpos + i];
        }
    }

    if (((MNStruct*)globalcontext)->storeFMatrices == 0) {

        for (int i = 0; i < FitRedCoeff / 2; i++) {
            for (int k = 0; k < ((MNStruct*)globalcontext)->pulse->nobs; k++) {
                double time = (double)((MNStruct*)globalcontext)->pulse->obsn[k].bat;

                TotalMatrix(k, i + TimetoMargin + startpos) = cos(2 * M_PI * freqs[i] * time);
                TotalMatrix(k, i + FitRedCoeff / 2 + TimetoMargin + startpos) =
                    sin(2 * M_PI * freqs[i] * time);
            }
        }
    }

    if (((MNStruct*)globalcontext)->incRED == 3) {

        for (int pl = 0; pl < ((MNStruct*)globalcontext)->numFitRedPL; pl++) {

            double Tspan = maxtspan;
            double f1yr = 1.0 / 3.16e7;

            double redamp = Cube[pcount];
            pcount++;
            double redindex = Cube[pcount];
            pcount++;

            double cornerfreq = 0;
            if (((MNStruct*)globalcontext)->incRED == 4) {
                cornerfreq = pow(10.0, Cube[pcount]) / Tspan;
                pcount++;
            }

            redamp = pow(10.0, redamp);
            if (((MNStruct*)globalcontext)->RedPriorType == 1) {
                uniformpriorterm += log(redamp);
            }

            double Agw = redamp;
            for (int i = 0; i < FitRedCoeff / 2; i++) {

                double rho = 0;
                if (((MNStruct*)globalcontext)->incRED == 3) {
                    rho = (Agw * Agw / 12.0 / (M_PI * M_PI)) * pow(f1yr, (-3)) *
                          pow(freqs[i] * 365.25, (-redindex)) / (Tspan * 24 * 60 * 60);
                }
                if (((MNStruct*)globalcontext)->incRED == 4) {

                    rho = pow((1 + (pow((1.0 / 365.25) / cornerfreq, redindex / 2))), 2) *
                          (Agw * Agw / 12.0 / (M_PI * M_PI)) /
                          pow((1 + (pow(freqs[i] / cornerfreq, redindex / 2))), 2) /
                          (Tspan * 24 * 60 * 60) * pow(f1yr, -3.0);
                }
                // if(rho > pow(10.0,15))rho=pow(10.0,15);
                powercoeff[i] += rho;
                powercoeff[i + FitRedCoeff / 2] += rho;
            }
        }

        startpos = FitRedCoeff;
    }

    for (int i = 0; i < FitRedCoeff / 2; i++) {
        freqdet = freqdet + 2 * log(powercoeff[i]);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)globalcontext)->incDM > 0) {

        for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {
            DMVec[o] = 1.0 / (DMKappa *
                              pow((double)((MNStruct*)globalcontext)->pulse->obsn[o].freqSSB, 2));
        }

        for (int i = 0; i < FitDMCoeff / 2; i++) {

            freqs[startpos + i] =
                ((MNStruct*)globalcontext)->sampleFreq[startpos / 2 + i] / maxtspan;
            freqs[startpos + i + FitDMCoeff / 2] = freqs[startpos + i];

            if (((MNStruct*)globalcontext)->storeFMatrices == 0) {
                for (int k = 0; k < ((MNStruct*)globalcontext)->pulse->nobs; k++) {
                    double time = (double)((MNStruct*)globalcontext)->pulse->obsn[k].bat;

                    TotalMatrix(k, (i + TimetoMargin + startpos)) =
                        cos(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
                    TotalMatrix(k, (i + FitDMCoeff / 2 + TimetoMargin + startpos)) =
                        sin(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
                }
            }
        }
    }

    if (((MNStruct*)globalcontext)->incDM == 3) {

        for (int pl = 0; pl < ((MNStruct*)globalcontext)->numFitDMPL; pl++) {
            double DMamp = Cube[pcount];
            pcount++;
            double DMindex = Cube[pcount];
            pcount++;

            double f1yr = 1.0 / 3.16e7;

            DMamp = pow(10.0, DMamp);
            if (((MNStruct*)globalcontext)->DMPriorType == 1) {
                uniformpriorterm += log(DMamp);
            }
            for (int i = 0; i < FitDMCoeff / 2; i++) {

                double rho = (DMamp * DMamp) * pow(f1yr, (-3)) *
                             pow(freqs[startpos + i] * 365.25, (-DMindex)) /
                             (maxtspan * 24 * 60 * 60);
                powercoeff[startpos + i] += rho;
                powercoeff[startpos + i + FitDMCoeff / 2] += rho;
            }
        }

        int coefftovary = 0;
        double amptovary = 0.0;

        for (int i = 0; i < FitDMCoeff / 2; i++) {
            freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
        }
        startpos += FitDMCoeff;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get Time domain likelihood//////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double tdet = 0;
    double timelike = 0;

    for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {
        timelike += Resvec[o] * Resvec[o] * Noise[o];
        tdet -= log(Noise[o]);
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////Do Algebra/////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    logtchk("Starting algebra");

    Eigen::MatrixXd NT = TotalMatrix;

    for (int i = 0; i < ((MNStruct*)globalcontext)->pulse->nobs; i++) {
        for (int j = 0; j < totalsize; j++) {
            NT(i, j) *= Noise[i];
        }
    }

    // std::cout << TotalMatrix << std::endl;

    Eigen::MatrixXd TNT = TotalMatrix.transpose() * NT;
    Eigen::VectorXd NTd = NT.transpose() * Resvec;

    logtchk("Finishing main algebra");

    for (int j = 0; j < totCoeff; j++) {
        TNT(TimetoMargin + j, (TimetoMargin + j)) += 1.0 / powercoeff[j];
    }

    // Perform Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(TNT);

    // Solve the linear system
    Eigen::VectorXd chol_solution = llt.solve(NTd);

    // Calculate log determinant
    double jointdet = 0;
    for (int i = 0; i < TNT.rows(); ++i) {
        jointdet += std::log(llt.matrixL()(i, i));
    }
    jointdet *= 2;  // Because det(A) = det(L)^2 for Cholesky A = LL^T

    double freqlike = 0;
    for (int j = 0; j < totalsize; j++) {
        freqlike += NTd[j] * chol_solution[j];
    }

    double lnewChol = -0.5 * (tdet + jointdet + freqdet + timelike - freqlike) + uniformpriorterm;

    // std::cout << "lnew " << lnewChol << " " << tdet << " " << jointdet << " " << freqdet << " "
    //          << timelike << " " << freqlike << " " << uniformpriorterm << std::endl;

    logtchk("Exiting TempoNest Likelihood");

    return lnewChol;
}
