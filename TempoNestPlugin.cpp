#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif
//  Copyright (C) 2013 Lindley Lentati

/*
 * This file is part of TempoNest
 *
 * TempoNest is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * TempoNest is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with TempoNest. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * If you use TempoNest and as a byproduct both Tempo2 and MultiNest
 * then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
 * Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2,
 * pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 * or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 * pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 * timing model and MultiNest Papers here.
 */

#include <dlfcn.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <t2fit.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include "TempoNest.h"
#include "multinest.h"
#include "tempo2.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"

#include <string.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

#include <gsl/gsl_sf_gamma.h>

#include <mpi.h>

#ifdef HAVE_MLAPACK
    #include <mpblas_mpfr.h>
#endif

void ephemeris_routines(pulsar* psr, int npsr);
void clock_corrections(pulsar* psr, int npsr);
void extra_delays(pulsar* psr, int npsr);

void fastephemeris_routines(pulsar* psr, int npsr)
{
    vectorPulsar(psr, npsr);     /* 1. Form a vector pointing at the pulsar */
    readEphemeris(psr, npsr, 0); /* 2. Read the ephemeris normally out */
    get_obsCoord(
        psr, npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre normally out*/
    tt2tb(psr, npsr); /* Observatory/time-dependent part of TT-TB normally out*/
    readEphemeris(psr, npsr, 0); /* Re-evaluate ephemeris with correct TB */
    // readShortEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */
}

void fastSubIntephemeris_routines(pulsar* psr, int npsr)
{
    // vectorPulsar(psr,npsr); /* 1. Form a vector pointing at the pulsar */
    // readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
    get_obsCoord(psr, npsr); /* 3. Get Coordinate of observatory relative to Earth's centre */
    // tt2tb(psr,npsr); /* Observatory/time-dependent part of TT-TB */
    // readEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */
}

void fastformBatsAll(pulsar* psr, int npsr)
{
    // clock_corrections(psr,npsr); /* Clock corrections ... */

    int dotime = 0;
    struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;
    if (dotime == 1) {
        gettimeofday(&tval_before, NULL);
    }

    fastephemeris_routines(psr, npsr); /* Ephemeris routines ... */

    if (dotime == 1) {

        gettimeofday(&tval_after, NULL);
        timersub(&tval_after, &tval_before, &tval_resultone);
        printf("Time elapsed Up to end of ephem: %ld.%06ld\n", (long int)tval_resultone.tv_sec,
               (long int)tval_resultone.tv_usec);
        gettimeofday(&tval_before, NULL);
    }

    extra_delays(psr, npsr); /* Other time delays ... */

    if (dotime == 1) {

        gettimeofday(&tval_after, NULL);
        timersub(&tval_after, &tval_before, &tval_resultone);
        printf("Time elapsed Up to end of extra: %ld.%06ld\n", (long int)tval_resultone.tv_sec,
               (long int)tval_resultone.tv_usec);
        gettimeofday(&tval_before, NULL);
    }

    formBats(psr, npsr); /* Form Barycentric arrival times */

    if (dotime == 1) {

        gettimeofday(&tval_after, NULL);
        timersub(&tval_after, &tval_before, &tval_resultone);
        printf("Time elapsed Up to end of formbats: %ld.%06ld\n", (long int)tval_resultone.tv_sec,
               (long int)tval_resultone.tv_usec);
        gettimeofday(&tval_before, NULL);
    }

    secularMotion(psr, npsr);
    // updateBatsAll(psr, npsr);
    if (dotime == 1) {

        gettimeofday(&tval_after, NULL);
        timersub(&tval_after, &tval_before, &tval_resultone);
        printf("Time elapsed Up to end of secular: %ld.%06ld\n", (long int)tval_resultone.tv_sec,
               (long int)tval_resultone.tv_usec);
        gettimeofday(&tval_before, NULL);
    }
}

void fastformSubIntBatsAll(pulsar* psr, int npsr)
{
    // clock_corrections(psr,npsr); /* Clock corrections ... */
    fastSubIntephemeris_routines(psr, npsr); /* Ephemeris routines ... */
    extra_delays(psr, npsr);                 /* Other time delays ... */
    formBats(psr, npsr);                     /* Form Barycentric arrival times */
    secularMotion(psr, npsr);
}

MNStruct* init_struct(pulsar* pulseval, long double** LDpriorsval, int numberpulsarsval,
                      int numFitJumpsval, int numFitTimingval, int systemcountval,
                      int numFitEFACval, int numFitEQUADval, int numFitRedCoeffval,
                      int numFitDMCoeffval, int numFitScatCoeffval, int numFitRedPLval,
                      int numFitDMPLval, int** TempoFitNumsval, int* TempoJumpNumsval,
                      int* sysFlagsval, int numdimsval, int incREDval, int incDMval,
                      int DMFloatstartval, int RedFloatstartval, int TimeMarginVal,
                      int JumpMarginVal, int doLinearVal, double* SampleFreqsVal,
                      char* whiteflagval, int whitemodelval, int RedPriorType, int DMPriorType,
                      int EQUADPriorType, int EFACPriorType, int useOriginalErrors,

                      int FitLowFreqCutoff, int useNbitsAlgebra, int sampler, int StoreFMatrices,
                      int debug, double* PreJumpVals, char* rootName, int doMax, int rank)
{
    MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));

    MNS->pulse = pulseval;
    MNS->LDpriors = LDpriorsval;
    MNS->numberpulsars = numberpulsarsval;
    MNS->numFitJumps = numFitJumpsval;
    MNS->numFitTiming = numFitTimingval;
    MNS->systemcount = systemcountval;
    MNS->numFitEFAC = numFitEFACval;
    MNS->numFitEQUAD = numFitEQUADval;
    MNS->numFitRedCoeff = numFitRedCoeffval;
    MNS->numFitDMCoeff = numFitDMCoeffval;
    MNS->numFitRedPL = numFitRedPLval;
    MNS->numFitDMPL = numFitDMPLval;
    MNS->TempoFitNums = TempoFitNumsval;
    MNS->TempoJumpNums = TempoJumpNumsval;
    MNS->sysFlags = sysFlagsval;
    MNS->numdims = numdimsval;
    MNS->incRED = incREDval;
    MNS->incDM = incDMval;

    MNS->TimeMargin = TimeMarginVal;
    MNS->JumpMargin = JumpMarginVal;
    MNS->doLinear = doLinearVal;
    MNS->sampleFreq = SampleFreqsVal;
    MNS->whiteflag = whiteflagval;
    MNS->whitemodel = whitemodelval;

    MNS->RedPriorType = RedPriorType;
    MNS->DMPriorType = DMPriorType;
    MNS->EQUADPriorType = EQUADPriorType;
    MNS->EFACPriorType = EFACPriorType;
    MNS->useOriginalErrors = useOriginalErrors;

    MNS->FitLowFreqCutoff = FitLowFreqCutoff;
    MNS->useNbitsAlgebra = useNbitsAlgebra;

    MNS->storeFMatrices = StoreFMatrices;

    MNS->debug = debug;
    MNS->PreJumpVals = PreJumpVals;
    MNS->rootName = rootName;
    MNS->doMax = doMax;
    MNS->Tspan = 0;
    MNS->TimetoMargin = 0;
    MNS->totCoeff = 0;
    MNS->totRedShapeCoeff = 0;
    MNS->totalsize = 0;

    MNS->rank = rank;
    return MNS;
}

void printPriors(pulsar* psr, long double** TempoPriors, double** Dpriors, int incEFAC,
                 int incEQUAD, int incRED, int incDM, int numRedCoeff, int numDMCoeff,
                 int numScatCoeff, int fitDMModel, std::string longname, void* context)
{

    std::ofstream getdistparamnames;
    std::string gdpnfname = longname + ".paramnames";
    std::string latexvar;
    getdistparamnames.open(gdpnfname.c_str());

    if (TempoPriors[0][2] == 0) {
        getdistparamnames << "Phase \\phi\n";
    }

    if (((MNStruct*)context)->rank == 0)
        printf("\nPriors:\n");
    int paramsfitted = 0;
    if (TempoPriors[0][2] == 0) {
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on Phase : %.25Lg -> %.25Lg\n",
                   TempoPriors[paramsfitted][0] +
                       Dpriors[paramsfitted][0] * TempoPriors[paramsfitted][1],
                   TempoPriors[paramsfitted][0] +
                       Dpriors[paramsfitted][1] * TempoPriors[paramsfitted][1]);
    }
    paramsfitted++;

    for (int p = 0; p < MAX_PARAMS; p++) {
        for (int k = 0; k < psr[0].param[p].aSize; k++) {
            if (p == param_raj)
                latexvar = "\\alpha";
            else if (p == param_decj)
                latexvar = "\\delta";
            else if (p == param_pmra)
                latexvar = "\\mu_\\alpha";
            else if (p == param_pmdec)
                latexvar = "\\mu_\\delta";
            else if (p == param_px)
                latexvar = "\\pi";
            else if (p == param_sini)
                latexvar = "\\sin i";
            else if (p == param_pb)
                latexvar = "P_b";
            else if (p == param_t0)
                latexvar = "T_0";
            else if (p == param_ecc)
                latexvar = "e";
            else if (p == param_pbdot)
                latexvar = "\\dot{P_bb}";
            else if (p == param_om)
                latexvar = "\\omega";
            else if (p == param_omdot)
                latexvar = "\\dot{\\omega}";
            else if (p == param_tasc)
                latexvar = "T_{asc}";
            else if (p == param_m2)
                latexvar = "m_2";
            else if (p == param_eps1)
                latexvar = "\\kappa";
            else if (p == param_eps2)
                latexvar = "\\eta";
            else if (p == param_gamma)
                latexvar = "\\gamma";
            else if (p == param_f && k == 0)
                latexvar = "\\nu";
            else if (p == param_f && k == 1)
                latexvar = "\\dot{\\nu}";
            else if (p == param_f && k == 2)
                latexvar = "\\ddot{\\nu}";
            else if (p == param_dm && k == 0)
                latexvar = "";
            else if (p == param_dm && k == 1)
                latexvar = "";
            else if (p == param_dm && k == 2)
                latexvar = "";
            else
                latexvar = psr[0].param[p].shortlabel[k];

            if (psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel) {
                if (TempoPriors[paramsfitted][2] == 0 && ((MNStruct*)context)->rank == 0) {
                    printf("Prior on %s : %.25Lg -> %.25Lg\n", psr[0].param[p].shortlabel[k],
                           TempoPriors[paramsfitted][0] +
                               Dpriors[paramsfitted][0] * TempoPriors[paramsfitted][1],
                           TempoPriors[paramsfitted][0] +
                               Dpriors[paramsfitted][1] * TempoPriors[paramsfitted][1]);
                }

                if (TempoPriors[paramsfitted][2] == 0) {
                    getdistparamnames << psr[0].param[p].shortlabel[k] << " " << latexvar << "\n";
                }
                paramsfitted++;
            }
        }
    }

    int jumpsfitted = 0;
    for (int i = 0; i <= psr[0].nJumps; i++) {
        if (psr[0].fitJump[i] == 1) {

            if (TempoPriors[paramsfitted][2] == 0) {

                if (((MNStruct*)context)->rank == 0)
                    printf("Prior on Jump %i : %.25Lg -> %.25Lg\n", jumpsfitted + 1,
                           TempoPriors[paramsfitted][0] +
                               Dpriors[paramsfitted][0] * TempoPriors[paramsfitted][1],
                           TempoPriors[paramsfitted][0] +
                               Dpriors[paramsfitted][1] * TempoPriors[paramsfitted][1]);

                getdistparamnames << "Jump" << jumpsfitted + 1 << "\n";
            }

            paramsfitted++;
            jumpsfitted++;
        }
    }

    if (incEFAC > 0) {
        int EFACnum = 1;
        for (int i = 0; i < incEFAC; i++) {
            if (((MNStruct*)context)->rank == 0)
                printf("Prior on EFAC %i : %.5g -> %.5g\n", EFACnum, Dpriors[paramsfitted][0],
                       Dpriors[paramsfitted][1]);
            getdistparamnames << "EFAC" << i + 1 << "\n";
            paramsfitted++;
            EFACnum++;
        }
    }

    if (incEQUAD > 0) {
        int EQUADnum = 1;
        for (int i = 0; i < incEQUAD; i++) {
            if (((MNStruct*)context)->rank == 0)
                printf("Prior on EQUAD %i: %.5g -> %.5g\n", EQUADnum, Dpriors[paramsfitted][0],
                       Dpriors[paramsfitted][1]);
            getdistparamnames << "EQUAD" << i + 1 << "\n";
            paramsfitted++;
            EQUADnum++;
        }
    }

    if (((MNStruct*)context)->FitLowFreqCutoff > 0) {
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on Red LF Cutoff %g -> %g \n", Dpriors[paramsfitted][0],
                   Dpriors[paramsfitted][1]);
        paramsfitted++;
        getdistparamnames << "LFC\n";
    }

    if (incRED == 1 || incRED == 3) {
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n", Dpriors[paramsfitted][0],
                   Dpriors[paramsfitted][1]);
        paramsfitted++;
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on Red Noise Slope : %.5g -> %.5g\n", Dpriors[paramsfitted][0],
                   Dpriors[paramsfitted][1]);
        paramsfitted++;

        getdistparamnames << "RedAmp\n";
        getdistparamnames << "RedSlope\n";

    } else if (incRED == 2) {
        int Coeffnum = 1;
        for (int i = 0; i < numRedCoeff; i++) {
            if (((MNStruct*)context)->rank == 0)
                printf("Prior on Red Noise Coefficient %i Log Amplitude : %.5g -> %.5g\n", Coeffnum,
                       Dpriors[paramsfitted][0], Dpriors[paramsfitted][1]);
            getdistparamnames << "RedC" << i + 1 << "\n";
            paramsfitted++;
        }
    }

    if (incDM == 1 || incDM == 3) {
        getdistparamnames << "DMAmp\n";
        getdistparamnames << "DMSlope\n";
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on DM Log Amplitude : %.5g -> %.5g\n", Dpriors[paramsfitted][0],
                   Dpriors[paramsfitted][1]);
        paramsfitted++;
        if (((MNStruct*)context)->rank == 0)
            printf("Prior on DM Slope : %.5g -> %.5g\n", Dpriors[paramsfitted][0],
                   Dpriors[paramsfitted][1]);
        paramsfitted++;

    } else if (incDM == 2) {
        int Coeffnum = 1;
        for (int i = 0; i < numDMCoeff; i++) {
            if (((MNStruct*)context)->rank == 0)
                printf("Prior on DM Coefficient %i Log Amplitude : %.5g -> %.5g\n", Coeffnum,
                       Dpriors[paramsfitted][0], Dpriors[paramsfitted][1]);
            getdistparamnames << "DMC" << i + 1 << "\n";
            paramsfitted++;
        }
    }

    if (fitDMModel == 1) {
        for (int i = 0; i < psr[0].dmoffsDMnum; i++) {
            getdistparamnames << "DMModel" << i + 1 << "\n";
            if (((MNStruct*)context)->rank == 0)
                printf("Prior for DMModel: %i %g %g \n", i, Dpriors[paramsfitted][0],
                       Dpriors[paramsfitted][1]);
            paramsfitted++;
        }
    }

    getdistparamnames.close();
}

void update_MNPriors(MNStruct* MNS, double** DPriorsval, long double** priorsval, int linearPriors)
{

    int paramsfitted = 1;

    if (linearPriors != 2) {
        // printf("lin p is 0\n");
        for (int p = 0; p < MAX_PARAMS; p++) {
            for (int k = 0; k < MNS->pulse->param[p].aSize; k++) {
                if (MNS->pulse->param[p].fitFlag[k] == 1 && p != param_dmmodel) {
                    if (p == param_ecc || p == param_px || p == param_m2 ||
                        p == param_dm && k == 0) {
                        long double minprior =
                            priorsval[paramsfitted][0] +
                            DPriorsval[paramsfitted + linearPriors][0] * priorsval[paramsfitted][1];
                        if (minprior < 0 && priorsval[paramsfitted][2] == 0) {
                            //				printf("%.10Lg %.10Lg %g %g
                            //\n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
                            long double newprior =
                                -priorsval[paramsfitted][0] / priorsval[paramsfitted][1];
                            DPriorsval[paramsfitted + linearPriors][0] = (double)newprior;
                            printf(
                                "Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n",
                                MNS->pulse->param[p].shortlabel[k],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][0] *
                                        priorsval[paramsfitted][1],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][1] *
                                        priorsval[paramsfitted][1]);
                            //              		printf("%.10Lg %.10Lg %g %g
                            //              \n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
                        }
                    }
                    if (p == param_stig) {
                        long double minprior =
                            priorsval[paramsfitted][0] +
                            DPriorsval[paramsfitted + linearPriors][0] * priorsval[paramsfitted][1];
                        long double maxprior =
                            priorsval[paramsfitted][0] +
                            DPriorsval[paramsfitted + linearPriors][1] * priorsval[paramsfitted][1];

                        if (minprior < 0 && priorsval[paramsfitted][2] == 0) {
                            long double newprior =
                                -priorsval[paramsfitted][0] / priorsval[paramsfitted][1];
                            DPriorsval[paramsfitted + linearPriors][0] = (double)newprior;
                            printf(
                                "Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n",
                                MNS->pulse->param[p].shortlabel[k],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][0] *
                                        priorsval[paramsfitted][1],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][1] *
                                        priorsval[paramsfitted][1]);
                        }
                        if (maxprior > 1 && priorsval[paramsfitted][2] == 0) {
                            long double newprior =
                                (1.0 - priorsval[paramsfitted][0]) / priorsval[paramsfitted][1];
                            DPriorsval[paramsfitted + linearPriors][1] = (double)newprior;
                            printf(
                                "Prior on %s updated to be physical (was >1) : %.25Lg -> %.25Lg\n",
                                MNS->pulse->param[p].shortlabel[k],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][0] *
                                        priorsval[paramsfitted][1],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][1] *
                                        priorsval[paramsfitted][1]);
                        }
                    }
                    if (p == param_sini) {

                        long double minprior =
                            priorsval[paramsfitted][0] +
                            DPriorsval[paramsfitted + linearPriors][0] * priorsval[paramsfitted][1];
                        long double maxprior =
                            priorsval[paramsfitted][0] +
                            DPriorsval[paramsfitted + linearPriors][1] * priorsval[paramsfitted][1];

                        if (minprior < 0 && priorsval[paramsfitted][2] == 0) {
                            long double newprior =
                                -priorsval[paramsfitted][0] / priorsval[paramsfitted][1];
                            DPriorsval[paramsfitted + linearPriors][0] = (double)newprior;
                            printf(
                                "Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n",
                                MNS->pulse->param[p].shortlabel[k],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][0] *
                                        priorsval[paramsfitted][1],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][1] *
                                        priorsval[paramsfitted][1]);
                        }
                        if (maxprior > 1 && priorsval[paramsfitted][2] == 0) {
                            long double newprior =
                                (1.0 - priorsval[paramsfitted][0]) / priorsval[paramsfitted][1];
                            DPriorsval[paramsfitted + linearPriors][1] = (double)newprior;
                            printf(
                                "Prior on %s updated to be physical (was >1) : %.25Lg -> %.25Lg\n",
                                MNS->pulse->param[p].shortlabel[k],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][0] *
                                        priorsval[paramsfitted][1],
                                priorsval[paramsfitted][0] +
                                    DPriorsval[paramsfitted + linearPriors][1] *
                                        priorsval[paramsfitted][1]);
                        }
                    }
                    paramsfitted++;
                }
            }
        }
    }

    MNS->Dpriors = DPriorsval;

    MNS->LDpriors = priorsval;
}

/************************************************* dumper routine
 * ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way
// he/she wants
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

void dumper(int& nSamples, int& nlive, int& nPar, double** physLive, double** posterior,
            double** paramConstr, double& maxLogLike, double& logZ, double& logZerr, void* context)
{
    // convert the 2D Fortran arrays to C++ arrays

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior
    // probability in the last two columns
    // printf("enter dump %i %i \n", nSamples, nPar);
    //	int i, j;
    //
    //	double **postdist=new double*[nSamples];
    //	for(int m1=0; m1<nSamples; m1++){
    //		postdist[m1]=new double[nPar+2];
    //	}
    //	for( i = 0; i < nPar + 2; i++ ){
    //		printf("%i\n",i);
    //		for( j = 0; j < nSamples; j++ )
    //			postdist[j][i] = posterior[0][i * nSamples + j];
    //	}

    // printf("mid dump \n");

    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last
    // column

    //	double pLivePts[nlive][nPar + 1];
    //	for( i = 0; i < nPar + 1; i++ )
    //		for( j = 0; j < nlive; j++ )
    //			pLivePts[j][i] = physLive[0][i * nlive + j];
    //
    ////	printf("exit dump\n");
    //	 for(int m1=0; m1<nSamples; m1++){
    //              free(postdist[m1]);
    //    }
    //	free(postdist);
}

/* The main function of a plugin called from Tempo2 is 'graphicalInterface'
 */
extern "C" int graphicalInterface(int argc, char** argv, pulsar* psr, int* pnpsr)
{
    int iteration;
    int listparms;
    int outRes = 0;
    int writeModel = 0;
    char timFile[MAX_PSR][MAX_FILELEN], parFile[MAX_PSR][MAX_FILELEN];
    char outputSO[MAX_FILELEN];
    char str[MAX_FILELEN];
    char newparname[MAX_FILELEN];
    longdouble coeff[MAX_COEFF]; /* For polynomial coefficients in polyco */
    int npsr = *pnpsr;           /* The number of pulsars */
    int noWarnings = 1;
    double globalParameter = 0.0;
    int displayParams;
    int nGlobal, i, flagPolyco = 0, it, k;
    char polyco_args[128];
    char polyco_file[128];
    int newpar = 0;
    int onlypre = 0;
    //  char tempo2MachineType[MAX_FILELEN]="";
    FILE* alias;
    char** commandLine;
    clock_t startClock, endClock;
    time_t rawstarttime, rawstoptime;
    struct tm* rawstarttimeinfo;
    struct tm* rawstoptimeinfo;
    const char* CVS_verNum = "$Revision: 1.28 $";
    int numFitJumps;
    int numToMargin = 0;
    int nprocs = 0;
    char* ConfigFileName;

    int rank, size;
    MPI_Comm world_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

    if (rank == 0) {
        printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
        printf("This is free software, and you are welcome to redistribute it\n");
        printf("under conditions of GPL license.\n\n");
    }

    startClock = clock();

    time(&rawstarttime);
    rawstarttimeinfo = localtime(&rawstarttime);
    time_t starttime = mktime(rawstarttimeinfo);
    if (rank == 0)
        printf("The start date/time is: %s", asctime(rawstarttimeinfo));
    commandLine = (char**)malloc(1000 * sizeof(char*));

    for (i = 0; i < 1000; i++)
        commandLine[i] = (char*)malloc(sizeof(char) * 1000);

    ConfigFileName = "defaultparameters.conf";
    /* Parse input line for machine type */
    for (i = 0; i < argc; i++) {
        if (strcasecmp(argv[i], "-Cfile") == 0) {
            ConfigFileName = argv[i + 1];
        }
        strcpy(commandLine[i], argv[i]);
    }

    if (rank == 0) {
        printf("Using config filename: %s \n", ConfigFileName);
    }

    strcpy(outputSO, "");
    npsr = 0; /* Initialise the number of pulsars */
    displayParams = 0;
    nGlobal = 0;
    /* Obtain command line arguments */
    logdbg("Running getInputs %d", psr[0].nits);
    logdbg("Completed getInputs");
    getInputs(psr, argc, commandLine, timFile, parFile, &listparms, &npsr, &nGlobal, &outRes,
              &writeModel, outputSO, &flagPolyco, polyco_args, polyco_file, &newpar, &onlypre,
              dcmFile, covarFuncFile, newparname);

    int doMaxLike = 0;

    logdbg("Reading par file");
    readParfile(psr, parFile, timFile,
                npsr); /* Read .par file to define the pulsar's initial parameters */
    logdbg("Finished reading par file %d", psr[0].nits);
    if (flagPolyco == 0) {
        logdbg("Running readTimfile");
        readTimfile(psr, timFile, npsr); /* Read .tim file to define the site-arrival-times */
        logdbg("Completed readTimfile %d", psr[0].param[param_ecc].paramSet[1]);
    }

    logdbg("Running preProcess %d", psr[0].nits);
    preProcess(psr, npsr, argc, commandLine);
    logdbg("Completed preProcess %d", psr[0].nits);

    if (debugFlag == 1) {
        logdbg("Number of iterations = %d", psr[0].nits);
        logdbg("Maximum number of parameters = %d", MAX_PARAMS);
        logdbg("Number of pulsars = %d", npsr);
    }

    int useGPUS;
    int useNbitsAlgebra = 0;

    int StoreFMatrices = 1;
    char root[100];
    int numTempo2its;
    int doLinearFit;
    int doMax;
    int incEFAC;
    int incEQUAD;

    int incRED;
    int FitLowFreqCutoff = 0;
    int incDM;

    int doTimeMargin;
    int doJumpMargin;
    double FitSig;
    int customPriors;
    int Reddims = 0;
    int DMdims = 0;
    int DMModeldims = 0;
    int fitDMModel = 0;
    double* EFACPrior;
    double* EQUADPrior;
    double* AlphaPrior;
    double* AmpPrior;
    double* DMAlphaPrior;
    double* DMAmpPrior;
    double numRedCoeff;
    double numDMCoeff;
    double numScatCoeff;
    int numRedPL;
    int numDMPL;
    double FourierSig;
    double* SampleFreq;
    int numEFAC = 0;
    int numEQUAD = 0;

    char wflag[100];
    int whitemodel;
    int RedPriorType;
    int DMPriorType;

    int EQUADPriorType;
    int EFACPriorType;
    int useOriginalErrors;

    char* Type = new char[100];
    char* WhiteName = new char[100];
    EFACPrior = new double[2];
    EQUADPrior = new double[2];
    AlphaPrior = new double[2];
    AmpPrior = new double[2];
    DMAlphaPrior = new double[2];
    DMAmpPrior = new double[2];

    int debug = 0;

    setupparams(ConfigFileName, useGPUS, Type, numTempo2its, doLinearFit, doMax, incEFAC, incEQUAD,
                incRED, incDM, doTimeMargin, doJumpMargin, FitSig, customPriors, EFACPrior,
                EQUADPrior, AlphaPrior, AmpPrior, DMAlphaPrior, DMAmpPrior, numRedCoeff, numDMCoeff,
                numScatCoeff, numRedPL, numDMPL, FourierSig, WhiteName, whitemodel, RedPriorType,
                DMPriorType, EQUADPriorType, EFACPriorType, useOriginalErrors, FitLowFreqCutoff,
                useNbitsAlgebra, StoreFMatrices, debug);

    if (rank == 0) {

#ifdef HAVE_MLAPACK
        if (useNbitsAlgebra > 0) {
            printf(
                "Using the MPLAPACK library with MPFR library for arbitrary decimal calculations");
            printf("MLAPACK: Set the linear algebra precision to %d bits\n", useNbitsAlgebra);
            mpfr_set_default_prec(useNbitsAlgebra);
        }
#endif
    }

    formBatsAll(psr, npsr); /* Form Barycentric arrival times */
    logdbg("calling formResiduals");
    formResiduals(psr, npsr, 1); /* Form residuals */

    /*Work out data time span to get maximum number of coefficients*/

    double start, end;
    int go = 0;
    for (int i = 0; i < psr[0].nobs; i++) {
        if (psr[0].obsn[i].deleted == 0) {
            if (go == 0) {
                go = 1;
                start = (double)psr[0].obsn[i].bat;
                end = start;
            } else {
                if (start > (double)psr[0].obsn[i].bat)
                    start = (double)psr[0].obsn[i].bat;
                if (end < (double)psr[0].obsn[i].bat)
                    end = (double)psr[0].obsn[i].bat;
            }
        }
    }

    double maxtspan = 1 * (end - start);

    if (maxtspan < 1) {
        maxtspan = maxtspan * 24 * 60;
        if (rank == 0)
            printf("Assume less than a day, Tspan is now in minutes\n");
    }

    int mindays = int(floor(1 + 2 * maxtspan / psr[0].nobs));
    int mincoeff = int(floor(1 + maxtspan / mindays));

    int Reddaysincoeffs = int(floor(maxtspan / numRedCoeff));
    int DMdaysincoeffs = int(floor(maxtspan / numDMCoeff));
    int Scatdaysincoeffs = int(floor(maxtspan / numScatCoeff));

    if (numRedCoeff < mindays) {
        numRedCoeff = mincoeff;
        if (FitLowFreqCutoff > 0) {
            numRedCoeff = mincoeff + 10;  // int(Reddaysincoeffs);
        }
    } else {
        numRedCoeff = int(Reddaysincoeffs);  // Reddaysincoeffs;
    }

    if (numDMCoeff < mindays) {
        numDMCoeff = mincoeff;
    } else {
        numDMCoeff = int(DMdaysincoeffs);  // DMdaysincoeffs;
    }

    if (numScatCoeff < mindays)
        numScatCoeff = mincoeff;
    else
        numScatCoeff = int(Scatdaysincoeffs);  // Scatdaysincoeffs;

    if (incRED == 0)
        numRedCoeff = 0;
    if (incDM == 0)
        numDMCoeff = 0;

    SampleFreq = new double[int(numRedCoeff + numDMCoeff + numScatCoeff)];
    setFrequencies(ConfigFileName, SampleFreq, numRedCoeff, numDMCoeff, numScatCoeff, 0, 0, 0, 1, 1,
                   1, 1, 1, 1);

    if (rank == 0)
        printf("Num T2 its %i \n", numTempo2its);
    for (it = 0; it < numTempo2its; it++) /* Why pulsar 0 should select the iterations? */
    {
        if (it > 0) /* Copy post-fit values to pre-fit values */
        {
            for (i = 0; i < MAX_PARAMS; i++) {
                for (int p = 0; p < npsr; p++) {
                    for (k = 0; k < psr[p].param[i].aSize; k++) {
                        psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
                        psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
                    }
                }
            }
        }
        //      long seed = TKsetSeed();
        for (iteration = 0; iteration < 2; iteration++) /* Do pre- and post- fit analysis */
        {
            logdbg("iteration %d", iteration);
            logdbg("calling formBatsAll");
            //	  printf("Calling formBats\n");
            formBatsAll(psr, npsr); /* Form Barycentric arrival times */
            logdbg("calling formResiduals");
            formResiduals(psr, npsr, 1); /* Form residuals */

            if (listparms == 1 && iteration == 0)
                displayParameters(13, timFile, parFile, psr,
                                  npsr); /* List out all the parameters */
            if (iteration == 0)          /* Only fit to pre-fit residuals */
            {
                logdbg("calling doFit");

                t2Fit(psr, npsr, covarFuncFile);
                // if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
                // doFit(psr,npsr,writeModel); /* Fit to the residuals to obtain updated parameters
                // */ else
                // doFitDCM(psr,dcmFile,covarFuncFile,npsr,writeModel);
                /* doFitGlobal(psr,npsr,&globalParameter,nGlobal,writeModel);*/ /* Fit to the
                                                                                   residuals to
                                                                                   obtain updated
                                                                                   parameters  */
                logdbg("completed doFit");
            }
            if (iteration == 1 || onlypre == 1) {
                if (strlen(outputSO) == 0) {
                    if (rank == 0)
                        textOutput(psr, npsr, globalParameter, nGlobal, outRes, newpar,
                                   newparname); /* Output results to the screen */
                } else                          /* Use a plug in for the output */
                {
                    char* (*entry)(int, char**, pulsar*, int);
                    void* module;
                    for (int iplug = 0; iplug < tempo2_plug_path_len; iplug++) {
                        sprintf(str, "%s/%s_%s_plug.t2", tempo2_plug_path[iplug], outputSO,
                                tempo2MachineType);
                        module = dlopen(str, RTLD_NOW);
                        if (module == NULL) {
                            printf("dlerror() = %s\n", dlerror());
                        } else
                            break;
                    }
                    if (!module) {
                        fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n");
                        return -1;
                    }
                    /*
                     * Check that the plugin is compiled against the same version of tempo2.h
                     */
                    char** pv = (char**)dlsym(module, "plugVersionCheck");
                    if (pv != NULL) {
                        // there is a version check for this plugin
                        if (strcmp(TEMPO2_h_VER, *pv)) {
                            fprintf(stderr, "[error]: Plugin version mismatch\n");
                            fprintf(stderr, " '%s' != '%s'\n", TEMPO2_h_VER, *pv);
                            fprintf(stderr,
                                    " Please recompile plugin against same tempo2 version!\n");
                            dlclose(module);
                            return -1;
                        }
                    }

                    entry = (char* (*)(int, char**, pulsar*, int))dlsym(module, "tempoOutput");
                    if (entry == NULL) {
                        dlclose(module);
                        fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n");
                        return -1;
                    }
                    entry(argc, argv, psr, npsr);
                }
            }
            psr[0].noWarnings = 2;
            if (onlypre == 1)
                iteration = 2;
            /*	  textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par");*/ /* Output
                                                                                             results
                                                                                             to the
                                                                                             screen
                                                                                           */
            /*	  printf("Next iteration\n");*/
        }
    }

    if (incRED == 0)
        Reddims = 0;
    if (incRED == 1)
        Reddims = 2;
    if (incRED == 2)
        Reddims = numRedCoeff;
    if (incRED == 3)
        Reddims = 2 * numRedPL;
    if (incRED == 4)
        Reddims = 3 * numRedPL;
    if (incRED == 5)
        Reddims = 2 * numRedCoeff + 2;

    if (FitLowFreqCutoff > 0)
        Reddims += 1;

    if (incDM == 0)
        DMdims = 0;
    if (incDM == 1)
        DMdims = 2;
    if (incDM == 2)
        DMdims = numDMCoeff;
    if (incDM == 3)
        DMdims = 2 * numDMPL;

    // printf("DMModel flag: %i \n",psr[0].param[param_dmmodel].fitFlag[0]);
    if (psr[0].param[param_dmmodel].fitFlag[0] == 1) {
        //	printf("Fitting for DMModel using %i structure functions \n",psr[0].dmoffsDMnum);
        DMModeldims = psr[0].dmoffsDMnum;
        fitDMModel = 1;
    }

    if (fitDMModel == 1 && (doTimeMargin != 0 || doJumpMargin != 0 || doLinearFit != 0)) {
        if (rank == 0)
            printf(
                "DMModel currently only supports doTimeMargin= 0, doJumpMargin = 0, doLinear = 0");
        // 			return 0;
    }

    if ((incRED == 1 && incDM != 1 && incDM != 0) || (incRED != 1 && incRED != 0 && incDM == 1)) {
        if (rank == 0)
            printf(
                "Different methods for DM and red noise not currently supported, please use the "
                "same option for both");
        return 0;
    }

    std::string pulsarname = psr[0].name;
    std::string longname = Type + pulsarname + "-";

    if (longname.size() >= 100) {
        if (rank == 0)
            printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n",
                   longname.size());
        return 0;
    }

    for (int r = 0; r <= longname.size(); r++) {
        root[r] = longname[r];
    }

    std::string wname = WhiteName;

    if (wname.size() >= 100) {
        if (rank == 0)
            printf(
                "white noise flag is too long, needs to be less than 100 characters, currently %i "
                ".\n",
                wname.size());
        return 0;
    }

    for (int r = 0; r <= wname.size(); r++) {
        wflag[r] = wname[r];
    }

    if (rank == 0) {
        printf("Graphical Interface: TempoNest\n");
        printf("Author:              L. Lentati\n");
        printf("Version:             1.0\n");
        printf("----------------------------------------------------------------\n");
        printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
        printf("This is free software, and you are welcome to redistribute it\n");
        printf("under conditions of GPL license.\n\n");
        printf("----------------------------------------------------------------\n");

        printf("\n\n\n\n*****************************************************\n");
        printf("Starting TempoNest\n");
        printf("*****************************************************\n\n\n\n");
        printf("Details of the fit:\n");
        printf("file root set to %s \n", root);
    }
    int systemcount = 0;
    int* numFlags = new int[psr[0].nobs];
    for (int o = 0; o < psr[0].nobs; o++) {
        numFlags[o] = 0;
    }

    for (int o = 0; o < psr[0].nobs; o++) {
        psr[0].obsn[o].snr = 1;
        psr[0].obsn[o].tobs = 1;
        for (int f = 0; f < psr[0].obsn[o].nFlags; f++) {
            if (strcasecmp(psr[0].obsn[o].flagID[f], "-snr") == 0) {
                psr[0].obsn[o].snr = atof(psr[0].obsn[o].flagVal[f]);
            }
            if (strcasecmp(psr[0].obsn[o].flagID[f], "-tobs") == 0) {
                psr[0].obsn[o].tobs = atof(psr[0].obsn[o].flagVal[f]);
            }
        }
    }

    //	if(incEFAC > 0 || incEQUAD > 0 || incShannonJitter > 0){printf("using white noise model
    //%i\n", whitemodel);}

    std::vector<std::string> systemnames;
    for (int o = 0; o < psr[0].nobs; o++) {
        int found = 0;
        for (int f = 0; f < psr[0].obsn[o].nFlags; f++) {
            if (strcasecmp(psr[0].obsn[o].flagID[f], wflag) == 0) {

                if (std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) !=
                    systemnames.end()) {
                    /* systemnames contains x */
                } else {
                    /* systemnames does not contain x */
                    if (incEFAC == 2 || incEQUAD == 2) {
                        if (rank == 0)
                            printf("Found %s %s \n", wflag, psr[0].obsn[o].flagVal[f]);
                    }
                    systemnames.push_back(psr[0].obsn[o].flagVal[f]);
                    systemcount++;
                }
                found = 1;
            }
        }
        if (found == 0 && (incEFAC == 2 || incEQUAD == 2)) {
            if (rank == 0)
                printf("Observation %i is missing the %s flag, please check before continuing\n", o,
                       wflag);
            return 0;
        }
    }

    if (incEFAC == 2 || incEQUAD == 2) {
        if (rank == 0)
            printf("total number of systems: %i \n", systemcount);
    }
    if (systemcount == 0 && (incEFAC < 2 && incEQUAD < 2)) {
        systemcount = 1;
    }

    for (int o = 0; o < psr[0].nobs; o++) {
        for (int f = 0; f < psr[0].obsn[o].nFlags; f++) {

            if (strcasecmp(psr[0].obsn[o].flagID[f], wflag) == 0) {
                for (int l = 0; l < systemcount; l++) {
                    if (psr[0].obsn[o].flagVal[f] == systemnames[l]) {
                        numFlags[o] = l;
                    }
                }
            }
        }
    }

    int* includeEQsys = new int[systemcount];
    for (int l = 0; l < systemcount; l++) {
        includeEQsys[l] = 1;
    }

    if (incEQUAD == 3) {
        for (int i = 0; i < psr[0].nTNEQ; i++) {
            if (rank == 0)
                printf("Including for EQUAD==3 %i %s \n", i, psr->TNEQFlagVal[i]);
            for (int l = 0; l < systemcount; l++) {
                if (psr->TNEQFlagVal[i] == systemnames[l]) {
                    if (rank == 0)
                        printf("This corresponds to system %i \n", l);
                    includeEQsys[l] = 0;
                }
            }
        }

        for (int l = 0; l < systemcount; l++) {
            int changed = 0;
            if (includeEQsys[l] == 0) {
                includeEQsys[l] = 1;
                changed = 1;
            }
            if (includeEQsys[l] == 1 && changed == 0) {
                includeEQsys[l] = 0;
            }
            if (rank == 0)
                printf("Including? %i %i \n", l, includeEQsys[l]);
        }
    }

    if (incEFAC == 1) {
        if (rank == 0)
            printf("Including One EFAC for all observations\n");
        numEFAC = 1;
    }
    if (incEFAC == 2) {
        if (rank == 0)
            printf("Including One EFAC for each %s\n", wflag);
        numEFAC = systemcount;
    }
    if (incEQUAD == 1) {
        if (rank == 0)
            printf("Including One EQUAD for all observations\n");
        numEQUAD = 1;
    }
    if (incEQUAD == 2) {
        if (rank == 0)
            printf("Including One EQUAD for each %s\n", wflag);
        numEQUAD = systemcount;
    }
    if (incEQUAD == 3) {
        if (rank == 0)
            printf("Including One EQUAD for the %i included TNEQ in par file\n", psr[0].nTNEQ);
        numEQUAD = psr[0].nTNEQ;
    }

    if (rank == 0) {

        if (incRED == 1)
            printf("Including Red Noise : Power Law Model\n");
        if (incRED == 2) {
            printf("Including Red Noise : Model Independant - Fitting %i Coefficients\n",
                   int(numRedCoeff));
        }
        if (incRED == 3) {
            printf("Including Red Noise: Power Law Model to %i Coefficients \n", int(numRedCoeff));
        }
        if (incRED == 4) {
            printf("Including Red Noise: Power Law Model with Corner Freq to %i Coefficients\n \n",
                   int(numRedCoeff));
        }
        if (incRED == 5) {
            printf("Including Red Noise Numerically: Power Law Model to %i Coefficients \n",
                   int(numRedCoeff));
        }
        if (incDM == 1) {
            printf("Including DM : Power Law Model\n");
        }
        if (incDM == 2) {
            printf("Including DM : Model Independant - Fitting %i Coefficients\n", int(numDMCoeff));
        }
        if (incDM == 3) {
            printf("Including DM: %i Component Power Law Model to %i Coefficients \n", numDMPL,
                   int(numDMCoeff));
        }
    }

    int fitcount = 0;
    if (rank == 0)
        printf("fitting for: Arbitrary Phase \n");
    fitcount++;
    for (int p = 0; p < MAX_PARAMS; p++) {
        for (int k = 0; k < psr[0].param[p].aSize; k++) {
            if (psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel) {
                if (rank == 0)
                    printf("fitting for: %s \n", psr[0].param[p].shortlabel[k]);
                fitcount++;
            }
            if (psr[0].param[p].fitFlag[k] == 1 && p == param_dmmodel) {
                if (rank == 0)
                    printf("fitting for: %s \n", psr[0].param[p].shortlabel[k]);
            }
        }
    }

    numFitJumps = 0;
    for (int i = 0; i <= psr[0].nJumps; i++) {
        if (psr[0].fitJump[i] == 1)
            numFitJumps++;
        // 		  printf("%i %i %g %g \n",i,psr[0].fitJump[i],
        // psr[0].jumpVal[i],psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree));
    }
    if (rank == 0) {
        printf("Found %i jumps to fit \n", numFitJumps);
        printf("total Timing Model params to fit:  %i \n", numFitJumps + fitcount);
    }

    long double** TempoPriors;
    long double* Tempo2Fit = new long double[numFitJumps + fitcount];
    int** TempoFitNums;
    int* TempoJumpNums;
    TempoPriors = new long double*[numFitJumps + fitcount];
    for (int i = 0; i < numFitJumps + fitcount; i++) {
        TempoPriors[i] = new long double[4];
        for (int j = 0; j < 4; j++) {
            TempoPriors[i][j] = 0;
        }
    }
    TempoFitNums = new int*[fitcount + numFitJumps];
    for (int i = 0; i < fitcount + numFitJumps; i++) {
        TempoFitNums[i] = new int[2];
        TempoFitNums[i][0] = -1;
        TempoFitNums[i][1] = -1;
    }
    TempoJumpNums = new int[numFitJumps];
    // printf("allocated\n");
    int paramsfitted = 0;
    //	printf("Offset and error: %g %g %g\n", psr[0].offset, psr[0].offset_e,
    // psr[0].offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree));
    TempoPriors[paramsfitted][0] = 0;  // psr[0].offset;
    TempoPriors[paramsfitted][1] =
        psr[0].offset_e;  // offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);

    /*
    if(GPTA == 1){
        TempoPriors[paramsfitted][0]=0;
            TempoPriors[paramsfitted][1]=1;//offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);
    }
    */

    TempoFitNums[paramsfitted][0] = 0;
    TempoFitNums[paramsfitted][1] = 0;
    // if(doTimeMargin != 0 || doJumpMargin != 0)TempoPriors[paramsfitted][2]=1;
    TempoPriors[paramsfitted][2] =
        1;  // always marginalise over phase. @TODO: make an option to allow fitting for
    paramsfitted++;
    for (int p = 0; p < MAX_PARAMS; p++) {
        for (int k = 0; k < psr[0].param[p].aSize; k++) {
            if (psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel) {

                TempoPriors[paramsfitted][0] = psr[0].param[p].prefit[k];
                TempoPriors[paramsfitted][1] =
                    psr[0].param[p].err[k] / std::sqrt(psr[0].fitChisq / psr[0].fitNfree);
                Tempo2Fit[paramsfitted] = psr[0].param[p].val[k];

                TempoFitNums[paramsfitted][0] = p;
                TempoFitNums[paramsfitted][1] = k;

                if (strcasecmp(psr[0].param[p].shortlabel[0], "F0") == 0) {
                    if (doTimeMargin != 0)
                        TempoPriors[paramsfitted][2] = 1;
                }
                if (doTimeMargin == 2)
                    TempoPriors[paramsfitted][2] = 1;

                paramsfitted++;
            }
        }
    }
    int jumpsfitted = 0;
    double* PreJumpVals = new double[psr[0].nJumps + 1]();
    if (rank == 0)
        printf("There are %i jumps \n", psr[0].nJumps);
    for (int i = 0; i <= psr[0].nJumps; i++) {

        char str1[100], str2[100], str3[100], str4[100], str5[100];
        int nread = sscanf(psr[0].jumpStr[i], "%s %s %s %s %s", str1, str2, str3, str4, str5);
        // printf("Jump String %s %s %s %s %s \n", str1, str2, str3, str4, str5);

        double prejump = 0;
        if (strcasecmp(str1, "MJD") == 0) {
            prejump = atof(str4);
            if (rank == 0)
                printf("MJD PreJump: %i %g \n", i, prejump);
        } else {
            prejump = atof(str3);
            if (rank == 0)
                printf("Not an MJD PreJump: %i %g \n", i, prejump);
        }
        PreJumpVals[i] = prejump;
        if (psr[0].fitJump[i] == 1) {

            //			printf("gonna read jump %i %s \n",i,psr[0].jumpStr[i]);
            char str1[100], str2[100], str3[100], str4[100], str5[100];
            int nread = sscanf(psr[0].jumpStr[i], "%s %s %s %s %s", str1, str2, str3, str4, str5);

            double prejump = 0;

            if (strcasecmp(str1, "MJD") == 0) {
                prejump = atof(str4);
            } else {
                prejump = atof(str3);
            }
            // printf(" Fitted Pre jump %i %g %g\n",i,prejump,psr[0].jumpVal[i]);
            // if(GPTA==0){
            if (numTempo2its > 1) {
                TempoPriors[paramsfitted][0] = psr[0].jumpVal[i];
                TempoPriors[paramsfitted][1] =
                    psr[0].jumpValErr[i] / std::sqrt(psr[0].fitChisq / psr[0].fitNfree);
                Tempo2Fit[paramsfitted] = psr[0].jumpVal[i];
            } else {
                TempoPriors[paramsfitted][0] = prejump;
                TempoPriors[paramsfitted][1] =
                    psr[0].jumpValErr[i] / std::sqrt(psr[0].fitChisq / psr[0].fitNfree);
                Tempo2Fit[paramsfitted] = prejump;
            }
            //}
            // else{

            //        TempoPriors[paramsfitted][0]=prejump;
            ///        TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
            //         Tempo2Fit[paramsfitted]=prejump;

            //}
            TempoJumpNums[jumpsfitted] = i;

            if (doJumpMargin == 1)
                TempoPriors[paramsfitted][2] = 1;

            paramsfitted++;
            jumpsfitted++;
        }
    }

    // set the MultiNest sampling parameters

    // 	return 0;
    int sampler = 1;
    int Nchords = 1;
    int IS = 1;        // do Nested Importance Sampling?
    int mmodal = 0;    // do mode separation?
    int ceff = 0;      // run in constant efficiency mode?
    int nlive = 500;   // number of live points
    double efr = 0.1;  // set the required efficiency
    int sample = 1;
    int nClsPar = 1;    // no. of parameters to do mode separation on
    int updInt = 2000;  // after how many iterations feedback is required & the output files should
                        // be updated
    int NBurn = 100;
    int NSamp = 10000;
    int GHSresume = 0;

    setupMNparams(ConfigFileName, sampler, IS, mmodal, ceff, nlive, efr, sample, updInt, nClsPar,
                  Nchords, NBurn, NSamp, GHSresume);

    double tol = 0.5;  // tol, defines the stopping criteria
    int ndims = numFitJumps + fitcount + numEFAC + numEQUAD + Reddims + DMdims +
                DMModeldims;  // dimensionality (no. of free parameters)
    int nPar = ndims;         // total no. of parameters including free & derived parameters
                       // note: posterior files are updated & dumper routine is called after every
                       // updInt*10 iterations
    double Ztol = -1E90;  // all the modes with logZ < Ztol are ignored
    int maxModes = 100;   // expected max no. of modes (used only for memory allocation)
    int pWrap[ndims];     // which parameters to have periodic boundary conditions?
    for (int i = 0; i < ndims; i++)
        pWrap[i] = 0;

    int seed = -1;    // random no. generator seed, if < 0 then take the seed from system clock
    int fb = 1;       // need feedback on standard output?
    int resume = 1;   // resume from a previous job?
    int outfile = 1;  // write output files?
    int initMPI = 0;  // initialize MPI routines?, relevant only if compiling with MPI set it to F
                      // if you want your main program to handle MPI initialization
    double logZero = -DBL_MAX;  // points with loglike < logZero will be ignored by MultiNest
    int maxiter = 0;  // max no. of iterations, a non-positive value means infinity. MultiNest will
                      // terminate if either it has done max no. of iterations or convergence
                      // criterion (defined through tol) has been satisfied
    void* context = 0;  // not required by MultiNest, any additional information user wants to pass
    // printf("Here \n");

    int FloatRedstart = numFitJumps + fitcount + numEFAC + numEQUAD + Reddims;
    int FloatDMstart = numFitJumps + fitcount + numEFAC + numEQUAD + Reddims + DMdims;

    char* chartroot = new char[longname.length() + 1];
    std::strcpy(chartroot, longname.c_str());

    MNStruct* MNS = init_struct(
        psr, TempoPriors, npsr, numFitJumps, fitcount, systemcount, numEFAC, numEQUAD,
        int(numRedCoeff), int(numDMCoeff), int(numScatCoeff), numRedPL, numDMPL, TempoFitNums,
        TempoJumpNums, numFlags, ndims, incRED, incDM, FloatDMstart, FloatRedstart, doTimeMargin,
        doJumpMargin, doLinearFit, SampleFreq, wflag, whitemodel, RedPriorType, DMPriorType,
        EQUADPriorType, EFACPriorType, useOriginalErrors, FitLowFreqCutoff, useNbitsAlgebra,
        sampler, StoreFMatrices, debug, PreJumpVals, chartroot, doMax, rank);

    MNS->includeEQsys = includeEQsys;

    // return 0;
    context = MNS;

    double** Dpriors;
    Dpriors = new double*[ndims];
    for (int i = 0; i < ndims; i++) {
        Dpriors[i] = new double[2];
    };

    // If using custompriors for errors incase T2 doesnt converge, get those values before doing
    // anything else
    if (customPriors == 1) {
        setTNPriors(ConfigFileName, Dpriors, TempoPriors,
                    ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps, ndims);
        update_MNPriors(MNS, Dpriors, TempoPriors, 0);

        context = MNS;

        int pcount = 1;
        for (int j = 1; j < ((MNStruct*)context)->numFitTiming; j++) {

            psr[0]
                .param[((MNStruct*)context)->TempoFitNums[j][0]]
                .val[((MNStruct*)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
            pcount++;
        }

        for (int j = 0; j < ((MNStruct*)context)->numFitJumps; j++) {
            psr[0].jumpVal[((MNStruct*)context)->TempoJumpNums[j]] = TempoPriors[pcount][0];
            pcount++;
        }
    }

    // Combine all the priors into one aray: Dpriors
    int pcount = 0;

    for (int i = 0; i < numFitJumps + fitcount; i++) {
        Dpriors[pcount][0] = -FitSig;
        Dpriors[pcount][1] = FitSig;
        pcount++;
    }

    for (int i = 0; i < numEFAC; i++) {
        Dpriors[pcount][0] = EFACPrior[0];
        Dpriors[pcount][1] = EFACPrior[1];
        pcount++;
    }
    for (int i = 0; i < numEQUAD; i++) {
        Dpriors[pcount][0] = EQUADPrior[0];
        Dpriors[pcount][1] = EQUADPrior[1];
        pcount++;
    }

    if (FitLowFreqCutoff > 0) {
        Dpriors[pcount][0] = -1;
        Dpriors[pcount][1] = 0;
        pcount++;
    }

    if (incRED == 1 || incRED == 3) {
        for (int i = 0; i < numRedPL; i++) {

            Dpriors[pcount][0] = AmpPrior[0];
            Dpriors[pcount][1] = AmpPrior[1];
            pcount++;
            Dpriors[pcount][0] = AlphaPrior[0];
            Dpriors[pcount][1] = AlphaPrior[1];
            pcount++;
        }
    }

    if (incDM == 1 || incDM == 3) {
        for (int i = 0; i < numDMPL; i++) {
            Dpriors[pcount][0] = DMAmpPrior[0];
            Dpriors[pcount][1] = DMAmpPrior[1];
            pcount++;
            Dpriors[pcount][0] = DMAlphaPrior[0];
            Dpriors[pcount][1] = DMAlphaPrior[1];
            pcount++;
        }
    }

    int* FitList = new int[ndims];
    for (int i = 0; i < ndims; i++) {
        FitList[i] = 0;
    }
    numToMargin = 0;
    for (int i = 0; i < ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps;
         i++) {
        if (TempoPriors[i][2] == 1) {
            FitList[i] = 1;
            if (rank == 0)
                printf("marginalising over timing param: %i \n", i);
            numToMargin++;
        }
    }

    getCustomDMatrix(psr, FitList, TempoFitNums, TempoJumpNums, Dpriors, incDM,
                     ((MNStruct*)context)->numFitTiming, ((MNStruct*)context)->numFitJumps);

    update_MNPriors(MNS, Dpriors, TempoPriors, 0);

    context = MNS;

    // Finally after doing everything reset custom priors in case overwritten by previous steps.
    if (customPriors == 1) {
        if (rank == 0)
            printf("Set to use custom priors, updating from setPriors function \n");
        setTNPriors(ConfigFileName, Dpriors, TempoPriors,
                    ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps, ndims);
        paramsfitted = 0;
    }
    update_MNPriors(MNS, Dpriors, TempoPriors, 0);
    context = MNS;

    printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff,
                numDMCoeff, numScatCoeff, fitDMModel, longname, context);

    if (rank == 0)
        printf("\n\n");
    ndims = ndims - numToMargin;
    nPar = ndims;

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////get TotalMatrix////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    getArraySizeInfo(context);
    if (rank == 0)
        printf("Pre-Computing Matrices: totalsize = %i\n", ((MNStruct*)context)->totalsize);

    for (int p = 1; p < ((MNStruct*)context)->numFitTiming; p++) {

        ((MNStruct*)context)
            ->pulse->param[((MNStruct*)context)->TempoFitNums[p][0]]
            .val[((MNStruct*)context)->TempoFitNums[p][1]] =
            ((MNStruct*)context)->Dpriors[p][0] * (((MNStruct*)context)->LDpriors[p][1]) +
            (((MNStruct*)context)->LDpriors[p][0]);
    }

    formBatsAll(((MNStruct*)context)->pulse, ((MNStruct*)context)->numberpulsars);
    formResiduals(((MNStruct*)context)->pulse, ((MNStruct*)context)->numberpulsars, 1);

    double* TotalMatrix =
        new double[((MNStruct*)context)->pulse->nobs * ((MNStruct*)context)->totalsize]();

    StoreTMatrix(TotalMatrix, context);

    ((MNStruct*)context)->StoredTMatrix = TotalMatrix;

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////call PolyChord/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double* PriorsArray = new double[2 * ndims];
    int NDerived = 0;

    int p = 0;
    pcount = 0;
    while (pcount < ndims) {
        if (((MNStruct*)context)->Dpriors[p][0] != ((MNStruct*)context)->Dpriors[p][1]) {
            PriorsArray[pcount] = ((MNStruct*)context)->Dpriors[p][0];
            PriorsArray[pcount + ndims] = ((MNStruct*)context)->Dpriors[p][1];
            if (rank == 0)
                printf("param %i in, %g %g \n", pcount, PriorsArray[pcount],
                       PriorsArray[pcount + ndims]);
            pcount++;
        }

        p++;
    }

    ((MNStruct*)context)->PriorsArray = PriorsArray;

    assigncontext(context);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////Call
    /// Samplers////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (sample == 1) {

        nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt,
                    Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter,
                    LRedLikeMNWrap, dumper, context);
    }

    if (rank == 0) {
        readsummary(psr, longname, ndims, context, Tempo2Fit, incRED, ndims, doTimeMargin,
                    doJumpMargin, doLinearFit);

        time(&rawstoptime);
        rawstoptimeinfo = localtime(&rawstoptime);
        time_t stoptime = mktime(rawstoptimeinfo);
        double seconds = difftime(stoptime, starttime);
        printf("The stop date/time was: %s", asctime(rawstoptimeinfo));

        printf("Total Wall clock run time: %g \n", seconds);
    }

    return EXIT_SUCCESS;
}

// redwards function to force linkage with library functions used by
// plugins
void thwart_annoying_dynamic_library_stuff(int never_call_me, float or_sink)
{
    ChebyModel* cm;
    T2Predictor* t2p;
    ChebyModel_Init(cm, 0, 0);
    T2Predictor_GetPhase(t2p, 0, 0);
}
