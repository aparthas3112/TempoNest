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

    double** EFAC;
    double* EQUAD;
    int pcount = 0;

    int numfit = ((MNStruct*)globalcontext)->numFitTiming + ((MNStruct*)globalcontext)->numFitJumps;
    int TimetoMargin = ((MNStruct*)globalcontext)->TimetoMargin;
    long double LDparams[numfit];
    for (int i = 0; i < numfit; i++) {
        LDparams[i] = 0;
    }
    double Fitparams[numfit];
    Eigen::VectorXd Resvec(((MNStruct*)globalcontext)->pulse->nobs);
    int fitcount = 0;

    pcount = 0;

    // Convert priors to physical units (here only for timing parameters and jumps)
    for (int p = 0;
         p < ((MNStruct*)globalcontext)->numFitTiming + ((MNStruct*)globalcontext)->numFitJumps;
         p++) {
        if (((MNStruct*)globalcontext)->Dpriors[p][1] !=
            ((MNStruct*)globalcontext)->Dpriors[p][0]) {
            double val = 0;
            if ((((MNStruct*)globalcontext)->LDpriors[p][3]) == 0) {
                val = Cube[fitcount];
            }
            if ((((MNStruct*)globalcontext)->LDpriors[p][3]) == 1) {
                val = pow(10.0, Cube[fitcount]);
            }
            if ((((MNStruct*)globalcontext)->LDpriors[p][3]) == 2) {
                val = pow(10.0, Cube[fitcount]);
                uniformpriorterm += log(val);
            }

            LDparams[p] = val * (((MNStruct*)globalcontext)->LDpriors[p][1]) +
                          (((MNStruct*)globalcontext)->LDpriors[p][0]);

            if (((MNStruct*)globalcontext)->TempoFitNums[p][0] == param_sini &&
                ((MNStruct*)globalcontext)->usecosiprior == 1) {
                val = Cube[fitcount];
                LDparams[p] = std::sqrt(1.0 - val * val);
            }

            fitcount++;

        } else if (((MNStruct*)globalcontext)->Dpriors[p][1] ==
                   ((MNStruct*)globalcontext)->Dpriors[p][0]) {
            LDparams[p] = ((MNStruct*)globalcontext)->Dpriors[p][0] *
                              (((MNStruct*)globalcontext)->LDpriors[p][1]) +
                          (((MNStruct*)globalcontext)->LDpriors[p][0]);
        }
    }

    pcount = 0;
    double phase = (double)LDparams[0];
    pcount++;
    for (int p = 1; p < ((MNStruct*)globalcontext)->numFitTiming; p++) {
        ((MNStruct*)globalcontext)
            ->pulse->param[((MNStruct*)globalcontext)->TempoFitNums[p][0]]
            .val[((MNStruct*)globalcontext)->TempoFitNums[p][1]] = LDparams[pcount];
        pcount++;
    }
    for (int p = 0; p < ((MNStruct*)globalcontext)->numFitJumps; p++) {
        ((MNStruct*)globalcontext)->pulse->jumpVal[((MNStruct*)globalcontext)->TempoJumpNums[p]] =
            LDparams[pcount];
        pcount++;
    }

    if (TimetoMargin != numfit) {
        fastformBatsAll(
            ((MNStruct*)globalcontext)->pulse,
            ((MNStruct*)globalcontext)->numberpulsars); /* Form Barycentric arrival times */
        formResiduals(((MNStruct*)globalcontext)->pulse, ((MNStruct*)globalcontext)->numberpulsars,
                      1); /* Form residuals */
    }

    for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {

        Resvec[o] = (double)((MNStruct*)globalcontext)->pulse->obsn[o].residual + phase;
    }

    pcount = fitcount;
    if (((MNStruct*)globalcontext)->incStep > 0) {
        for (int i = 0; i < ((MNStruct*)globalcontext)->incStep; i++) {

            int GrouptoFit = 0;
            double GroupStartTime = 0;

            GrouptoFit = floor(Cube[pcount]);
            pcount++;

            double GLength = ((MNStruct*)globalcontext)->GroupStartTimes[GrouptoFit][1] -
                             ((MNStruct*)globalcontext)->GroupStartTimes[GrouptoFit][0];
            GroupStartTime =
                ((MNStruct*)globalcontext)->GroupStartTimes[GrouptoFit][0] + Cube[pcount] * GLength;
            pcount++;

            double StepAmp = Cube[pcount];
            pcount++;

            // printf("Step details: Group %i S %g F %g SS %g A %g \n", GrouptoFit, ((MNStruct
            // *)globalcontext)->GroupStartTimes[GrouptoFit][0],((MNStruct
            // *)globalcontext)->GroupStartTimes[GrouptoFit][1],GroupStartTime,StepAmp);

            for (int o1 = 0; o1 < ((MNStruct*)globalcontext)->pulse->nobs; o1++) {
                if (((MNStruct*)globalcontext)->pulse->obsn[o1].sat > GroupStartTime &&
                    ((MNStruct*)globalcontext)->GroupNoiseFlags[o1] == GrouptoFit) {
                    Resvec[o1] += StepAmp;
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Subtract GLitches/////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < ((MNStruct*)globalcontext)->incGlitch; i++) {
        double GlitchMJD = Cube[pcount];
        pcount++;

        double* GlitchAmps = new double[3];
        if (((MNStruct*)globalcontext)->incGlitchTerms == 1) {
            GlitchAmps[0] = Cube[pcount];
            pcount++;
        } else if (((MNStruct*)globalcontext)->incGlitchTerms == 2) {
            GlitchAmps[0] = Cube[pcount];
            pcount++;
            GlitchAmps[1] = Cube[pcount];
            pcount++;
        } else if (((MNStruct*)globalcontext)->incGlitchTerms == 3) {
            GlitchAmps[0] = Cube[pcount];
            pcount++;
            GlitchAmps[1] = pow(10.0, Cube[pcount]);  // Decay Amp
            pcount++;
            GlitchAmps[2] = pow(10.0, Cube[pcount]);  // Decay Timescale
            pcount++;
        }

        for (int o1 = 0; o1 < ((MNStruct*)globalcontext)->pulse->nobs; o1++) {
            if (((MNStruct*)globalcontext)->pulse->obsn[o1].bat > GlitchMJD) {

                if (((MNStruct*)globalcontext)->incGlitchTerms == 1) {

                    long double arg = 0;
                    arg = ((((MNStruct*)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
                           ((MNStruct*)globalcontext)->pulse->param[param_f].val[0]) *
                          86400.0;
                    double darg = (double)arg;
                    Resvec[o1] += GlitchAmps[0] * darg;

                } else if (((MNStruct*)globalcontext)->incGlitchTerms == 2) {
                    for (int j = 0; j < ((MNStruct*)globalcontext)->incGlitchTerms; j++) {

                        long double arg = 0;
                        if (j == 0) {
                            arg = ((((MNStruct*)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
                                   ((MNStruct*)globalcontext)->pulse->param[param_f].val[0]) *
                                  86400.0;
                        }
                        if (j == 1) {
                            arg =
                                0.5 *
                                pow((((MNStruct*)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) *
                                        86400.0,
                                    2) /
                                ((MNStruct*)globalcontext)->pulse->param[param_f].val[0];
                        }
                        double darg = (double)arg;
                        Resvec[o1] += GlitchAmps[j] * darg;
                    }
                } else if (((MNStruct*)globalcontext)->incGlitchTerms == 3) {
                    long double arg = 0;
                    double time = (((MNStruct*)globalcontext)->pulse->obsn[o1].bat - GlitchMJD);
                    arg = ((((MNStruct*)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
                           ((MNStruct*)globalcontext)->pulse->param[param_f].val[0]) *
                          86400.0;
                    double darg = (double)arg;
                    Resvec[o1] += GlitchAmps[0] * darg +
                                  GlitchAmps[1] * darg * exp(-1 * time / GlitchAmps[2]);
                }
            }
        }

        delete[] GlitchAmps;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Get White Noise vector///////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)globalcontext)->numFitEFAC == 0) {
        EFAC = new double*[((MNStruct*)globalcontext)->EPolTerms];
        for (int n = 1; n <= ((MNStruct*)globalcontext)->EPolTerms; n++) {
            EFAC[n - 1] = new double[((MNStruct*)globalcontext)->systemcount];
            if (n == 1) {
                for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
                    EFAC[n - 1][o] = 1;
                }
            } else {
                for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
                    EFAC[n - 1][o] = 0;
                }
            }
        }
    } else if (((MNStruct*)globalcontext)->numFitEFAC == 1) {
        EFAC = new double*[((MNStruct*)globalcontext)->EPolTerms];
        for (int n = 1; n <= ((MNStruct*)globalcontext)->EPolTerms; n++) {

            EFAC[n - 1] = new double[((MNStruct*)globalcontext)->systemcount];
            if (n == 1) {
                for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
                    EFAC[n - 1][o] = pow(10.0, Cube[pcount]);
                    if (((MNStruct*)globalcontext)->EFACPriorType == 1) {
                        uniformpriorterm += log(EFAC[n - 1][o]);
                    }
                }
                pcount++;
            } else {
                for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {

                    EFAC[n - 1][o] = pow(10.0, Cube[pcount]);
                }
                pcount++;
            }
        }

    } else if (((MNStruct*)globalcontext)->numFitEFAC > 1) {
        EFAC = new double*[((MNStruct*)globalcontext)->EPolTerms];
        for (int n = 1; n <= ((MNStruct*)globalcontext)->EPolTerms; n++) {
            EFAC[n - 1] = new double[((MNStruct*)globalcontext)->systemcount];
            if (n == 1) {
                for (int p = 0; p < ((MNStruct*)globalcontext)->systemcount; p++) {
                    EFAC[n - 1][p] = pow(10.0, Cube[pcount]);
                    if (((MNStruct*)globalcontext)->EFACPriorType == 1) {
                        uniformpriorterm += log(EFAC[n - 1][p]);
                    }
                    pcount++;
                }
            } else {
                for (int p = 0; p < ((MNStruct*)globalcontext)->systemcount; p++) {
                    EFAC[n - 1][p] = pow(10.0, Cube[pcount]);
                    pcount++;
                }
            }
        }
    }

    // printf("Equad %i \n", ((MNStruct *)globalcontext)->numFitEQUAD);
    if (((MNStruct*)globalcontext)->numFitEQUAD == 0) {
        EQUAD = new double[((MNStruct*)globalcontext)->systemcount];
        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
            EQUAD[o] = 0;
        }
    } else if (((MNStruct*)globalcontext)->numFitEQUAD == 1) {
        EQUAD = new double[((MNStruct*)globalcontext)->systemcount];
        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {
            EQUAD[o] = pow(10.0, 2 * Cube[pcount]);
            if (((MNStruct*)globalcontext)->EQUADPriorType == 1) {
                uniformpriorterm += log(pow(10.0, Cube[pcount]));
            }
        }
        pcount++;
    } else if (((MNStruct*)globalcontext)->numFitEQUAD > 1) {
        EQUAD = new double[((MNStruct*)globalcontext)->systemcount];
        for (int o = 0; o < ((MNStruct*)globalcontext)->systemcount; o++) {

            if (((MNStruct*)globalcontext)->includeEQsys[o] == 1) {
                // printf("Cube: %i %i %g \n", o, pcount, Cube[pcount]);
                EQUAD[o] = pow(10.0, 2 * Cube[pcount]);
                if (((MNStruct*)globalcontext)->EQUADPriorType == 1) {
                    uniformpriorterm += log(pow(10.0, Cube[pcount]));
                }
                pcount++;
            } else {
                EQUAD[o] = 0;
            }
            // printf("Equad? %i %g \n", o, EQUAD[o]);
        }
    }

    double* Noise;
    double* BATvec;
    Noise = new double[((MNStruct*)globalcontext)->pulse->nobs];
    // BATvec=new double[((MNStruct *)globalcontext)->pulse->nobs];

    // for(int o=0;o<((MNStruct *)globalcontext)->pulse->nobs; o++){
    // BATvec[o]=(double)((MNStruct *)globalcontext)->pulse->obsn[o].bat;
    //}

    double DMKappa = 2.410 * pow(10.0, -16);
    if (((MNStruct*)globalcontext)->whitemodel == 0) {

        for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {
            double EFACterm = 0;
            double noiseval = 0;

            if (((MNStruct*)globalcontext)->useOriginalErrors == 0) {
                noiseval = ((MNStruct*)globalcontext)->pulse->obsn[o].toaErr;
            } else if (((MNStruct*)globalcontext)->useOriginalErrors == 1) {
                noiseval = ((MNStruct*)globalcontext)->pulse->obsn[o].origErr;
            }

            for (int n = 1; n <= ((MNStruct*)globalcontext)->EPolTerms; n++) {
                EFACterm =
                    EFACterm + pow((noiseval * pow(10.0, -6)) / pow(pow(10.0, -7), n - 1), n) *
                                   EFAC[n - 1][((MNStruct*)globalcontext)->sysFlags[o]];
            }

            // printf("Noise: %i %g %g %g %g %g \n", EFACterm, EQUAD[((MNStruct
            // *)globalcontext)->sysFlags[o]], ShannonJitterTerm, SWTerm, DMEQUADTerm);
            Noise[o] = 1.0 / (pow(EFACterm, 2) + EQUAD[((MNStruct*)globalcontext)->sysFlags[o]]);
        }

    } else if (((MNStruct*)globalcontext)->whitemodel == 1) {

        for (int o = 0; o < ((MNStruct*)globalcontext)->pulse->nobs; o++) {

            Noise[o] =
                1.0 /
                (EFAC[0][((MNStruct*)globalcontext)->sysFlags[o]] *
                 EFAC[0][((MNStruct*)globalcontext)->sysFlags[o]] *
                 (pow(((((MNStruct*)globalcontext)->pulse->obsn[o].toaErr) * pow(10.0, -6)), 2) +
                  EQUAD[((MNStruct*)globalcontext)->sysFlags[o]]));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////Initialise TotalMatrix////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    int totalsize = ((MNStruct*)globalcontext)->totalsize;

    Eigen::MatrixXd TotalMatrix(((MNStruct*)context)->pulse->nobs, ((MNStruct*)context)->totalsize);

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double maxtspan = ((MNStruct*)globalcontext)->Tspan;
    double averageTSamp = 2 * maxtspan / ((MNStruct*)globalcontext)->pulse->nobs;

    int FitRedCoeff = 2 * (((MNStruct*)globalcontext)->numFitRedCoeff);
    int FitDMCoeff = 2 * (((MNStruct*)globalcontext)->numFitDMCoeff);

    int totCoeff = ((MNStruct*)globalcontext)->totCoeff;

    double* powercoeff = new double[totCoeff];
    for (int o = 0; o < totCoeff; o++) {
        powercoeff[o] = 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double* freqs = new double[totCoeff];
    double* DMVec = new double[((MNStruct*)globalcontext)->pulse->nobs];

    double freqdet = 0;
    int startpos = 0;

    if (((MNStruct*)globalcontext)->incRED > 0 || ((MNStruct*)globalcontext)->incGWB == 1) {

        if (((MNStruct*)globalcontext)->FitLowFreqCutoff == 1) {
            double fLow = pow(10.0, Cube[pcount]);
            pcount++;

            double deltaLogF = 0.1;
            double RedMidFreq = 2.0;

            double RedLogDiff = log10(RedMidFreq) - log10(fLow);
            int LogLowFreqs = floor(RedLogDiff / deltaLogF);

            double RedLogSampledDiff = LogLowFreqs * deltaLogF;
            double sampledFLow = floor(log10(fLow) / deltaLogF) * deltaLogF;

            int freqStartpoint = 0;

            for (int i = 0; i < LogLowFreqs; i++) {
                ((MNStruct*)globalcontext)->sampleFreq[freqStartpoint] =
                    pow(10.0, sampledFLow + i * RedLogSampledDiff / LogLowFreqs);
                freqStartpoint++;
            }

            for (int i = 0; i < FitRedCoeff / 2 - LogLowFreqs; i++) {
                ((MNStruct*)globalcontext)->sampleFreq[freqStartpoint] = i + RedMidFreq;
                freqStartpoint++;
            }
        }

        if (((MNStruct*)globalcontext)->FitLowFreqCutoff == 2) {
            double fLow = pow(10.0, Cube[pcount]);
            pcount++;

            for (int i = 0; i < FitRedCoeff / 2; i++) {
                ((MNStruct*)globalcontext)->sampleFreq[i] = ((double)(i + 1)) * fLow;
            }
        }

        for (int i = 0; i < FitRedCoeff / 2 - ((MNStruct*)globalcontext)->incFloatRed; i++) {

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

    if (((MNStruct*)globalcontext)->incRED == 2) {

        for (int i = 0; i < FitRedCoeff / 2; i++) {
            int pnum = pcount;
            double pc = Cube[pcount];

            if (((MNStruct*)globalcontext)->RedPriorType == 1) {
                uniformpriorterm += log(pow(10.0, pc));
            }

            powercoeff[i] = pow(10.0, 2 * pc);
            powercoeff[i + FitRedCoeff / 2] = powercoeff[i];
            pcount++;
        }

        startpos = FitRedCoeff;

    } else if (((MNStruct*)globalcontext)->incRED == 3 || ((MNStruct*)globalcontext)->incRED == 4) {

        for (int pl = 0; pl < ((MNStruct*)globalcontext)->numFitRedPL; pl++) {

            double Tspan = maxtspan;
            double f1yr = 1.0 / 3.16e7;

            if (((MNStruct*)globalcontext)->FitLowFreqCutoff == 2) {
                Tspan = Tspan / ((MNStruct*)globalcontext)->sampleFreq[0];
            }

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
            for (int i = 0; i < FitRedCoeff / 2 - ((MNStruct*)globalcontext)->incFloatRed; i++) {

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

        int coefftovary = 0;
        double amptovary = 0.0;
        if (((MNStruct*)globalcontext)->varyRedCoeff == 1) {
            coefftovary = int(pow(10.0, Cube[pcount])) - 1;
            pcount++;
            amptovary = pow(10.0, Cube[pcount]);
            pcount++;

            powercoeff[coefftovary] = amptovary;
            powercoeff[coefftovary + FitRedCoeff / 2] = amptovary;
        }

        startpos = FitRedCoeff;
    }

    if (((MNStruct*)globalcontext)->incGWB == 1) {
        double GWBAmp = pow(10.0, Cube[pcount]);
        pcount++;
        uniformpriorterm += log(GWBAmp);
        double Tspan = maxtspan;

        if (((MNStruct*)globalcontext)->FitLowFreqCutoff == 2) {
            Tspan = Tspan / ((MNStruct*)globalcontext)->sampleFreq[0];
        }

        double f1yr = 1.0 / 3.16e7;
        for (int i = 0; i < FitRedCoeff / 2 - ((MNStruct*)globalcontext)->incFloatRed; i++) {
            double rho = (GWBAmp * GWBAmp / 12.0 / (M_PI * M_PI)) * pow(f1yr, (-3)) *
                         pow(freqs[i] * 365.25, (-4.333)) / (Tspan * 24 * 60 * 60);
            powercoeff[i] += rho;
            powercoeff[i + FitRedCoeff / 2] += rho;
            // printf("%i %g %g \n", i, freqs[i], powercoeff[i]);
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
                ((MNStruct*)globalcontext)
                    ->sampleFreq[startpos / 2 - ((MNStruct*)globalcontext)->incFloatRed + i] /
                maxtspan;
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

    if (((MNStruct*)globalcontext)->incDM == 2) {

        for (int i = 0; i < FitDMCoeff / 2; i++) {
            int pnum = pcount;
            double pc = Cube[pcount];

            powercoeff[startpos + i] = pow(10.0, 2 * pc);
            powercoeff[startpos + i + FitDMCoeff / 2] = powercoeff[startpos + i];
            freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
            pcount++;
        }
        startpos += FitDMCoeff;

    } else if (((MNStruct*)globalcontext)->incDM == 3) {

        for (int pl = 0; pl < ((MNStruct*)globalcontext)->numFitDMPL; pl++) {
            double DMamp = Cube[pcount];
            pcount++;
            double DMindex = Cube[pcount];
            pcount++;

            double Tspan = maxtspan;
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
        if (((MNStruct*)globalcontext)->varyDMCoeff == 1) {
            coefftovary = int(pow(10.0, Cube[pcount])) - 1;
            pcount++;
            amptovary = pow(10.0, Cube[pcount]) / (maxtspan * 24 * 60 * 60);
            pcount++;

            powercoeff[startpos + coefftovary] = amptovary;
            powercoeff[startpos + coefftovary + FitDMCoeff / 2] = amptovary;
        }

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
        // printf("Res: %i  %g %g \n", o, Resvec[o], sqrt(Noise[o]));
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

    delete[] DMVec;

    for (int i = 0; i < ((MNStruct*)globalcontext)->EPolTerms; i++)
        delete[] EFAC[i];
    delete[] EFAC;
    delete[] EQUAD;
    delete[] powercoeff;
    delete[] freqs;
    delete[] Noise;

    ((MNStruct*)globalcontext)->PreviousJointDet = jointdet;
    ((MNStruct*)globalcontext)->PreviousFreqDet = freqdet;
    ((MNStruct*)globalcontext)->PreviousUniformPrior = uniformpriorterm;

    // printf("tdet %g, jointdet %g, freqdet %g, lnew %g, timelike %g, freqlike %g\n", tdet,
    // jointdet, freqdet, lnew, timelike, freqlike);

    // printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
    logtchk("Exiting TempoNest Likelihood");

    return lnewChol;
}
