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
#include "eigen_config.h"
#include "namespaces/sampler.h"
#include "namespaces/settings.h"
#include "types/model.h"

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
    struct timeval tval_before, tval_after, tval_resultone;
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

MNStruct* init_struct(pulsar* pulseval, int numberpulsarsval, int timing_model_params,
                      int systemcountval, int numFitRedCoeffval, int numFitDMCoeffval,
                      int* sysFlagsval, int numdimsval, char* whiteflagval, int useOriginalErrors,
                      int debug, char* rootName, int rank)
{
    MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));

    MNS->pulse = pulseval;
    MNS->numberpulsars = numberpulsarsval;
    MNS->TimetoMargin = timing_model_params;

    MNS->systemcount = systemcountval;
    MNS->numFitRedCoeff = numFitRedCoeffval;
    MNS->numFitDMCoeff = numFitDMCoeffval;

    MNS->sysFlags = sysFlagsval;
    MNS->numdims = numdimsval;

    MNS->whiteflag = whiteflagval;

    MNS->useOriginalErrors = useOriginalErrors;

    MNS->debug = debug;
    MNS->rootName = rootName;
    MNS->Tspan = 0;
    MNS->totCoeff = 0;
    MNS->totalsize = 0;

    MNS->rank = rank;
    return MNS;
}

void printPriors(std::string longname)
{

    std::ofstream getdistparamnames;
    std::string gdpnfname = longname + ".paramnames";
    getdistparamnames.open(gdpnfname.c_str());

    if (model::efac.has_value()) {

        getdistparamnames << "EFAC\n";
    }

    if (model::equad.has_value()) {

        getdistparamnames << "EQUAD\n";
    }

    if (model::pl_red_noise) {

        getdistparamnames << "RedAmp\n";
        getdistparamnames << "RedSlope\n";
    }

    if (model::pl_dm_noise.has_value()) {
        getdistparamnames << "DMAmp\n";
        getdistparamnames << "DMSlope\n";
    }

    getdistparamnames.close();
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
    int npsr = *pnpsr; /* The number of pulsars */
    double globalParameter = 0.0;
    int nGlobal, i, flagPolyco = 0, it, k;
    char polyco_args[128];
    char polyco_file[128];
    int newpar = 0;
    int onlypre = 0;
    char** commandLine;
    time_t rawstarttime, rawstoptime;
    struct tm* rawstarttimeinfo;
    struct tm* rawstoptimeinfo;
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
    nGlobal = 0;
    /* Obtain command line arguments */
    logdbg("Running getInputs %d", psr[0].nits);
    logdbg("Completed getInputs");
    getInputs(psr, argc, commandLine, timFile, parFile, &listparms, &npsr, &nGlobal, &outRes,
              &writeModel, outputSO, &flagPolyco, polyco_args, polyco_file, &newpar, &onlypre,
              dcmFile, covarFuncFile, newparname);

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

    settings::load_settings("model_config.json");
    model::load_model("model_config.json");

    char root[100];
    int numTempo2its;
    int incEFAC;
    int incEQUAD;

    int Reddims = 0;
    int DMdims = 0;
    double numRedCoeff;
    double numDMCoeff;

    char wflag[100];

    int useOriginalErrors;

    char* Type = new char[100];
    char* WhiteName = new char[100];

    int debug = 0;

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

    if (numRedCoeff < mindays) {
        numRedCoeff = mincoeff;

    } else {
        numRedCoeff = int(Reddaysincoeffs);  // Reddaysincoeffs;
    }

    if (numDMCoeff < mindays) {
        numDMCoeff = mincoeff;
    } else {
        numDMCoeff = int(DMdaysincoeffs);  // DMdaysincoeffs;
    }

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
        }
    }

    std::string pulsarname = psr[0].name;
    std::string longname = Type + pulsarname + "-";

    if (longname.size() >= 100) {
        if (rank == 0)
            printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n",
                   (int)longname.size());
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
                (int)wname.size());
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

    int timing_model_params = 0;
    timing_model_params++;
    for (int p = 0; p < MAX_PARAMS; p++) {
        for (int k = 0; k < psr[0].param[p].aSize; k++) {
            if (psr[0].param[p].fitFlag[k] == 1) {
                if (rank == 0)
                    printf("fitting for: %s \n", psr[0].param[p].shortlabel[k]);
                timing_model_params++;
            }
        }
    }

    for (int i = 0; i <= psr[0].nJumps; i++) {
        if (psr[0].fitJump[i] == 1)
            timing_model_params++;
    }

    // set the MultiNest sampling parameters

    sampler::load_sampler("model_config.json");

    double tol = 0.5;  // tol, defines the stopping criteria
    int ndims = 4;

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

    char* chartroot = new char[longname.length() + 1];
    std::strcpy(chartroot, longname.c_str());

    MNStruct* MNS =
        init_struct(psr, npsr, timing_model_params, systemcount, int(numRedCoeff), int(numDMCoeff),
                    numFlags, ndims, wflag, useOriginalErrors, debug, chartroot, rank);

    context = MNS;

    printPriors(longname);

    if (rank == 0)
        printf("\n\n");

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////get TotalMatrix////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    getArraySizeInfo(context);
    if (rank == 0)
        printf("Pre-Computing Matrices: totalsize = %i\n", ((MNStruct*)context)->totalsize);

    formBatsAll(((MNStruct*)context)->pulse, ((MNStruct*)context)->numberpulsars);
    formResiduals(((MNStruct*)context)->pulse, ((MNStruct*)context)->numberpulsars, 1);

    double* TotalMatrix =
        new double[((MNStruct*)context)->pulse->nobs * ((MNStruct*)context)->totalsize]();

    StoreTMatrix(TotalMatrix, context);

    ((MNStruct*)context)->StoredTMatrix = TotalMatrix;

    std::cout << "Stored t matrix " << ndims << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////call PolyChord/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    assigncontext(context);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////Call
    /// Samplers////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (sampler::sample) {

        std::cout << "run " << std::endl;

        nested::run(sampler::importance_sampling, sampler::modal, sampler::constant_efficiency,
                    sampler::live_points, tol, sampler::efficiency, ndims, ndims,
                    sampler::num_cluster_parameters, maxModes, sampler::update_interval, Ztol, root,
                    seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLikeMNWrap,
                    dumper, context);
    }

    if (rank == 0) {
        readsummary(psr, longname, ndims, context, ndims);

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
