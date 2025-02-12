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
#include "samplers/multinest_interface.h"
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
#include "pulsar_utils.h"
#include "tests/tests.h"
#include "types/model.h"

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

void dumper(int& nSamples, int& nlive, int& nPar, double** physLive, double** posterior, double** paramConstr, double& maxLogLike, double& logZ, double& logZerr, void* context) {}

/* The main function of a plugin called from Tempo2 is 'graphicalInterface'
 */
extern "C" int graphicalInterface(int argc, char** argv, pulsar* psr, int* pnum_pulsars)
{
    int iteration;
    int listparms;
    int outRes = 0;
    int writeModel = 0;
    char timFile[MAX_PSR][MAX_FILELEN], parFile[MAX_PSR][MAX_FILELEN];
    char outputSO[MAX_FILELEN];
    char str[MAX_FILELEN];
    char newparname[MAX_FILELEN];
    int num_pulsars = *pnum_pulsars; /* The number of pulsars */
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

    ConfigFileName = "defaultparameters.json";
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
    num_pulsars = 0; /* Initialise the number of pulsars */
    nGlobal = 0;
    /* Obtain command line arguments */
    logdbg("Running getInputs %d", globals::pulsar->nits);
    logdbg("Completed getInputs");
    getInputs(psr, argc, commandLine, timFile, parFile, &listparms, &num_pulsars, &nGlobal, &outRes, &writeModel, outputSO, &flagPolyco, polyco_args, polyco_file, &newpar, &onlypre, dcmFile,
              covarFuncFile, newparname);

    logdbg("Reading par file");
    readParfile(psr, parFile, timFile, num_pulsars); /* Read .par file to define the pulsar's initial parameters */
    logdbg("Finished reading par file %d", globals::pulsar->nits);
    if (flagPolyco == 0) {
        logdbg("Running readTimfile");
        readTimfile(psr, timFile, num_pulsars); /* Read .tim file to define the site-arrival-times */
        logdbg("Completed readTimfile %d", globals::pulsar->param[param_ecc].paramSet[1]);
    }

    std::cout << "call pre process" << std::endl;
    logdbg("Running preProcess %d", globals::pulsar->nits);
    preProcess(psr, num_pulsars, argc, commandLine);
    logdbg("Completed preProcess %d", globals::pulsar->nits);

    std::cout << "initialise pulsar" << std::endl;

    globals::pulsar = &psr[0];
    initialise_pulsar(onlypre);

    std::cout << "load settings" << std::endl;

    globals::load_settings(ConfigFileName);
    model::load_model(ConfigFileName);
    sampler::load_sampler(ConfigFileName);

    std::string pulsarname = globals::pulsar->name;
    std::string longname = globals::root + pulsarname + "-";

    if (longname.size() >= 100) {
        if (rank == 0)
            printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n", (int)longname.size());
        return 0;
    }

    char root[100];
    for (int r = 0; r <= longname.size(); r++) {
        root[r] = longname[r];
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

    for (int o = 0; o < globals::pulsar->nobs; o++) {
        globals::pulsar->obsn[o].snr = 1;
        globals::pulsar->obsn[o].tobs = 1;
        for (int f = 0; f < globals::pulsar->obsn[o].nFlags; f++) {
            if (strcasecmp(globals::pulsar->obsn[o].flagID[f], "-snr") == 0) {
                globals::pulsar->obsn[o].snr = atof(globals::pulsar->obsn[o].flagVal[f]);
            }
            if (strcasecmp(globals::pulsar->obsn[o].flagID[f], "-tobs") == 0) {
                globals::pulsar->obsn[o].tobs = atof(globals::pulsar->obsn[o].flagVal[f]);
            }
        }
    }

    // set the MultiNest sampling parameters

    double tol = 0.5;  // tol, defines the stopping criteria
    int ndims = model::model_space.get_fitted_dims();

    double Ztol = -1E90;  // all the modes with logZ < Ztol are ignored
    int maxModes = 100;   // expected max no. of modes (used only for memory allocation)
    int pWrap[ndims];     // which parameters to have periodic boundary conditions?
    for (int i = 0; i < ndims; i++)
        pWrap[i] = 0;

    int seed = -1;              // random no. generator seed, if < 0 then take the seed from system clock
    int fb = 1;                 // need feedback on standard output?
    int resume = 1;             // resume from a previous job?
    int outfile = 1;            // write output files?
    int initMPI = 0;            // initialize MPI routines?, relevant only if compiling with MPI set it to F
                                // if you want your main program to handle MPI initialization
    double logZero = -DBL_MAX;  // points with loglike < logZero will be ignored by MultiNest
    int maxiter = 0;            // max no. of iterations, a non-positive value means infinity. MultiNest will
                                // terminate if either it has done max no. of iterations or convergence
                                // criterion (defined through tol) has been satisfied
    void* context = 0;          // not required by MultiNest, any additional information user wants to pass
    // printf("Here \n");

    char* chartroot = new char[longname.length() + 1];
    std::strcpy(chartroot, longname.c_str());

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////get TotalMatrix////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    model::model_space.update_array_size_info();

    formBatsAll(globals::pulsar, num_pulsars);
    formResiduals(globals::pulsar, num_pulsars, 1);

    model::model_space.store_total_matrix();

    // if we are running unit tests do that now rather than sampling
    if (globals::test_mode) {
        run_tests();
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////Call
    /// Samplers////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (sampler::sample) {

        std::cout << "run " << std::endl;

        nested::run(sampler::importance_sampling, sampler::modal, sampler::constant_efficiency, sampler::live_points, tol, sampler::efficiency, ndims, ndims, sampler::num_cluster_parameters, maxModes,
                    sampler::update_interval, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLikeMNWrap, dumper, 0);
    }

    if (rank == 0) {
        readsummary(globals::pulsar, longname, ndims, 0, ndims);

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
