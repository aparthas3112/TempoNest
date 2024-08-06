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

#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "T2toolkit.h"
#include "TempoNest.h"
#include "tempo2.h"

int UtWrap(int kX, int const kLowerBound, int const kUpperBound)
{
    int range_size = kUpperBound - kLowerBound + 1;

    if (kX < kLowerBound)
        kX += range_size * ((kLowerBound - kX) / range_size + 1);

    return kLowerBound + (kX - kLowerBound) % range_size;
}

double iter_factorial(unsigned int n)
{
    double ret = 1;
    for (unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

void readtxtoutput(std::string longname, int ndim, double** paramarray)
{

    int number_of_lines = 0;
    double weightsum = 0;

    std::ifstream checkfile;
    std::string checkname = longname + ".txt";
    checkfile.open(checkname.c_str());
    std::string line;
    while (getline(checkfile, line))
        ++number_of_lines;

    //        printf("number of lines %i \n",number_of_lines);
    checkfile.close();

    std::ifstream summaryfile;
    std::string fname = longname + ".txt";
    summaryfile.open(fname.c_str());

    for (int i = 0; i < ndim; i++) {
        paramarray[i][0] = 0;
        paramarray[i][1] = 0;
    }

    printf("Processing Posteriors \n");
    //	printf("Getting Means \n");
    double maxlike = -1.0 * pow(10.0, 10);
    double MAP = 0;
    for (int i = 0; i < number_of_lines; i++) {

        std::string line;
        getline(summaryfile, line);
        std::istringstream myStream(line);
        std::istream_iterator<double> begin(myStream), eof;
        std::vector<double> paramlist(begin, eof);

        weightsum += paramlist[0];
        double like = paramlist[1];

        if (like > maxlike) {
            maxlike = like;
            for (int i = 0; i < ndim; i++) {
                paramarray[i][2] = paramlist[i + 2];
            }
        }

        if (paramlist[0] > MAP) {
            MAP = paramlist[0];
            for (int i = 0; i < ndim; i++) {
                paramarray[i][3] = paramlist[i + 2];
            }
        }

        for (int i = 0; i < ndim; i++) {
            paramarray[i][0] += paramlist[i + 2] * paramlist[0];
        }
    }

    for (int i = 0; i < ndim; i++) {
        paramarray[i][0] = paramarray[i][0] / weightsum;
    }

    summaryfile.close();

    summaryfile.open(fname.c_str());

    printf("Getting Errors \n");

    for (int i = 0; i < number_of_lines; i++) {

        std::string line;
        getline(summaryfile, line);
        std::istringstream myStream(line);
        std::istream_iterator<double> begin(myStream), eof;
        std::vector<double> paramlist(begin, eof);

        for (int i = 0; i < ndim; i++) {
            paramarray[i][1] += paramlist[0] * (paramlist[i + 2] - paramarray[i][0]) *
                                (paramlist[i + 2] - paramarray[i][0]);
        }
    }

    for (int i = 0; i < ndim; i++) {
        paramarray[i][1] = sqrt(paramarray[i][1] / weightsum);
    }

    summaryfile.close();
}

void readphyslive(std::string longname, int ndim, double** paramarray, int sampler)
{

    int number_of_lines = 0;

    std::ifstream checkfile;
    std::string checkname;
    if (sampler == 0) {
        checkname = longname + "phys_live.points";
    }
    if (sampler == 1) {
        checkname = longname + "_phys_live.txt";
    }
    checkfile.open(checkname.c_str());
    std::string line;
    while (getline(checkfile, line))
        ++number_of_lines;

    checkfile.close();

    std::ifstream summaryfile;
    std::string fname;
    if (sampler == 0) {
        fname = longname + "phys_live.points";
    }
    if (sampler == 1) {
        fname = longname + "_phys_live.txt";
    }

    summaryfile.open(fname.c_str());

    double* TempCube = new double[ndim];
    for (int i = 0; i < ndim; i++) {
        paramarray[i][2] = 0;
    }

    //	printf("Getting ML \n");
    double maxlike = -1.0 * pow(10.0, 10);
    for (int i = 0; i < number_of_lines; i++) {

        std::string line;
        getline(summaryfile, line);
        std::istringstream myStream(line);
        std::istream_iterator<double> begin(myStream), eof;
        std::vector<double> paramlist(begin, eof);

        double like = paramlist[ndim];

        if (like > maxlike) {
            maxlike = like;
            for (int i = 0; i < ndim; i++) {
                paramarray[i][2] = paramlist[i];
                TempCube[i] = paramlist[i];
            }
        }
    }
    summaryfile.close();

    // void *tempcontext;
    // double *DerivedParams=new double[1];
    // int np = 1;
    // double like = NewLRedMarginLogLike(ndim, TempCube, np, DerivedParams, tempcontext);
}

void read_ddgr(std::string longname, int ndim, int sampler, void* context)
{

    int number_of_lines = 0;

    std::ifstream checkfile;
    std::string checkname;
    if (sampler == 0) {
        checkname = longname + "post_equal_weights.dat";
    }
    if (sampler == 1) {
        checkname = longname + "_phys_live.txt";
    }
    checkfile.open(checkname.c_str());
    std::string line;
    while (getline(checkfile, line))
        ++number_of_lines;

    checkfile.close();

    std::ifstream summaryfile;
    std::string fname;
    if (sampler == 0) {
        fname = longname + "post_equal_weights.dat";
    }
    if (sampler == 1) {
        fname = longname + "_phys_live.txt";
    }

    summaryfile.open(fname.c_str());

    double* TempCube = new double[ndim];

    //	printf("Getting ML \n");
    double maxlike = -1.0 * pow(10.0, 10);
    for (int i = 0; i < number_of_lines; i++) {

        std::string line;
        getline(summaryfile, line);
        std::istringstream myStream(line);
        std::istream_iterator<double> begin(myStream), eof;
        std::vector<double> paramlist(begin, eof);

        double like = paramlist[ndim];

        for (int i = 0; i < ndim; i++) {
            TempCube[i] = paramlist[i];
        }
        long double
            LDparams[((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps];
        int fitcount = 0;
        for (int p = 0; p < ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps;
             p++) {
            if (((MNStruct*)context)->Dpriors[p][1] != ((MNStruct*)context)->Dpriors[p][0]) {

                double val = 0;
                if ((((MNStruct*)context)->LDpriors[p][3]) == 0) {
                    val = TempCube[fitcount];
                }
                if ((((MNStruct*)context)->LDpriors[p][3]) == 1) {
                    val = pow(10.0, TempCube[fitcount]);
                }
                if ((((MNStruct*)context)->LDpriors[p][3]) == 2) {
                    val = pow(10.0, TempCube[fitcount]);
                }
                LDparams[p] = val * (((MNStruct*)context)->LDpriors[p][1]) +
                              (((MNStruct*)context)->LDpriors[p][0]);

                if (((MNStruct*)context)->TempoFitNums[p][0] == param_sini &&
                    ((MNStruct*)context)->usecosiprior == 1) {
                    val = TempCube[fitcount];
                    LDparams[p] = sqrt(1.0 - val * val);
                }

                fitcount++;

            } else if (((MNStruct*)context)->Dpriors[p][1] == ((MNStruct*)context)->Dpriors[p][0]) {
                LDparams[p] =
                    ((MNStruct*)context)->Dpriors[p][0] * (((MNStruct*)context)->LDpriors[p][1]) +
                    (((MNStruct*)context)->LDpriors[p][0]);
            }
        }
        int pcount = 0;
        pcount++;
        for (int p = 1; p < ((MNStruct*)context)->numFitTiming; p++) {
            ((MNStruct*)context)
                ->pulse->param[((MNStruct*)context)->TempoFitNums[p][0]]
                .val[((MNStruct*)context)->TempoFitNums[p][1]] = LDparams[pcount];
            pcount++;
        }
        for (int p = 0; p < ((MNStruct*)context)->numFitJumps; p++) {
            ((MNStruct*)context)->pulse->jumpVal[((MNStruct*)context)->TempoJumpNums[p]] =
                LDparams[pcount];
            pcount++;
        }

        fastformBatsAll(((MNStruct*)context)->pulse,
                        ((MNStruct*)context)->numberpulsars); /* Form Barycentric arrival times */
        formResiduals(((MNStruct*)context)->pulse, ((MNStruct*)context)->numberpulsars,
                      1); /* Form residuals */
        DDGRmodel(((MNStruct*)context)->pulse, 0, 0, -2);

        printf("PBDot %.20Lg  OMDot %.20Lg Gamma %.20Lg \n",
               ((MNStruct*)context)->pulse->param[param_pbdot].val[0],
               ((MNStruct*)context)->pulse->param[param_omdot].val[0],
               ((MNStruct*)context)->pulse->param[param_gamma].val[0]);
    }
    summaryfile.close();

    // void *tempcontext;
    // double *DerivedParams=new double[1];
    // int np = 1;
    // double like = NewLRedMarginLogLike(ndim, TempCube, np, DerivedParams, tempcontext);
}

void readsummary(pulsar* psr, std::string longname, int ndim, void* context, long double* Tempo2Fit,
                 int incRED, int ndims, int doTimeMargin, int doJumpMargin, int doLinear)
{

    int pcount = 0;
    int fitcount = 0;

    std::vector<double> paramlist(2 * ndims);

    double** paramarray = new double*[ndims];
    for (int p = 0; p < ndims; p++) {
        paramarray[p] = new double[4];
    }

    readtxtoutput(longname, ndim, paramarray);
    readphyslive(longname, ndim, paramarray, ((MNStruct*)context)->sampler);
    //	read_ddgr( longname, ndim, ((MNStruct *)context)->sampler,  context);
    int numlongparams = ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps;
    long double* LDP = new long double[numlongparams];
    pcount = 0;
    fitcount = 0;
    for (int j = 0; j < ((MNStruct*)context)->numFitTiming; j++) {
        if (((MNStruct*)context)->Dpriors[pcount][0] != ((MNStruct*)context)->Dpriors[pcount][1]) {
            double val = 0;
            if (((MNStruct*)context)->LDpriors[pcount][3] == 0) {
                val = paramarray[fitcount][2];
            }
            if (((MNStruct*)context)->LDpriors[pcount][3] == 1) {
                val = pow(10.0, paramarray[fitcount][2]);
            }

            LDP[j] = val * (((MNStruct*)context)->LDpriors[pcount][1]) +
                     (((MNStruct*)context)->LDpriors[pcount][0]);

            if (((MNStruct*)context)->TempoFitNums[pcount][0] == param_sini &&
                ((MNStruct*)context)->usecosiprior == 1) {
                val = paramarray[fitcount][2];
                LDP[j] = sqrt(1.0 - val * val);
            }

            fitcount++;
        } else if (((MNStruct*)context)->Dpriors[pcount][0] ==
                   ((MNStruct*)context)->Dpriors[pcount][1]) {
            LDP[j] = ((MNStruct*)context)->Dpriors[pcount][0] *
                         (((MNStruct*)context)->LDpriors[pcount][1]) +
                     (((MNStruct*)context)->LDpriors[pcount][0]);
        }
        //	printf("LD: %.25Lg %g %g \n",LDP[j], ((MNStruct
        //*)context)->Dpriors[pcount][0],((MNStruct *)context)->Dpriors[pcount][1]);
        pcount++;
    }

    for (int j = 0; j < ((MNStruct*)context)->numFitJumps; j++) {
        if (((MNStruct*)context)->Dpriors[pcount][0] != ((MNStruct*)context)->Dpriors[pcount][1]) {
            LDP[pcount] = paramarray[fitcount][2] * (((MNStruct*)context)->LDpriors[pcount][1]) +
                          (((MNStruct*)context)->LDpriors[pcount][0]);
            fitcount++;
        } else if (((MNStruct*)context)->Dpriors[pcount][0] ==
                   ((MNStruct*)context)->Dpriors[pcount][1]) {
            LDP[pcount] = ((MNStruct*)context)->Dpriors[pcount][0] *
                              (((MNStruct*)context)->LDpriors[pcount][1]) +
                          (((MNStruct*)context)->LDpriors[pcount][0]);
        }
        pcount++;
    }

    pcount = 0;
    int jj = 0;
    pcount = 1;
    fitcount = 0;
    if (((MNStruct*)context)->LDpriors[0][2] == 0)
        fitcount++;
    for (int j = 1; j < ((MNStruct*)context)->numFitTiming; j++) {

        long double value;
        long double error;

        value = LDP[pcount];

        if (((MNStruct*)context)->LDpriors[j][2] == 0) {
            error = paramarray[fitcount][1] * (((MNStruct*)context)->LDpriors[pcount][1]);
            fitcount++;
        } else if (((MNStruct*)context)->LDpriors[j][2] == 1) {
            error = 0;
        }

        ((MNStruct*)context)
            ->pulse->param[((MNStruct*)context)->TempoFitNums[pcount][0]]
            .val[((MNStruct*)context)->TempoFitNums[pcount][1]] = value;
        ((MNStruct*)context)
            ->pulse->param[((MNStruct*)context)->TempoFitNums[pcount][0]]
            .err[((MNStruct*)context)->TempoFitNums[pcount][1]] = error;
        pcount++;
    }

    for (int j = 0; j < ((MNStruct*)context)->numFitJumps; j++) {

        long double value;
        long double error;

        value = LDP[pcount];

        if (((MNStruct*)context)->LDpriors[pcount][2] == 0) {
            error = paramarray[fitcount][1] * (((MNStruct*)context)->LDpriors[pcount][1]);
            fitcount++;
        } else if (((MNStruct*)context)->LDpriors[pcount][2] == 1) {
            error = 0;
        }

        ((MNStruct*)context)->pulse->jumpVal[((MNStruct*)context)->TempoJumpNums[j]] = value;
        ((MNStruct*)context)->pulse->jumpValErr[((MNStruct*)context)->TempoJumpNums[j]] = error;
        pcount++;
    }

    if (((MNStruct*)context)->incBreakingIndex == 1) {
        long double F0 = ((MNStruct*)context)->pulse->param[param_f].val[0];
        long double F1 = ((MNStruct*)context)->pulse->param[param_f].val[1];
        long double BIndex = (long double)paramarray[fitcount][2];
        fitcount++;
        long double Kappa = -F1 / pow(F0, BIndex);  //(long double)pow(10.0, Cube[pcount]);
        fitcount++;
        long double F2 = -Kappa * BIndex * F1 * pow(F0, BIndex - 1);

        printf("%.15Lg %.15Lg %.15Lg \n", F0, F1, F2);
        ((MNStruct*)context)->pulse->param[param_f].val[2] = F2;
    }

    formBatsAll(((MNStruct*)context)->pulse, 1);  // Form Barycentric arrival times
    // printf("formed bats \n");
    formResiduals(((MNStruct*)context)->pulse, 1, 1);  // Form residuals
    // printf("done bats and stuff \n");
    std::ofstream designfile;
    std::string dname = longname + "T2scaling.txt";

    designfile.open(dname.c_str());
    double pdParamDeriv[MAX_PARAMS];
    int numtofit = ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps;
    for (int i = 1; i < ((MNStruct*)context)->numFitTiming; i++) {
        designfile << psr->param[((MNStruct*)context)->TempoFitNums[i][0]]
                          .label[((MNStruct*)context)->TempoFitNums[i][1]];
        designfile << " ";
        std::stringstream ss;
        ss.precision(std::numeric_limits<long double>::digits);  // override the default

        ss << ((MNStruct*)context)->LDpriors[i][0];
        ss << " ";
        ss << ((MNStruct*)context)->LDpriors[i][1];
        designfile << ss.str();
        designfile << "\n";
    }

    designfile.close();
    printf("text output\n");
    double Evidence = 0;
    TNtextOutput(((MNStruct*)context)->pulse, 1, 0, Tempo2Fit, context, incRED, ndims, paramlist,
                 Evidence, doTimeMargin, doJumpMargin, doLinear, longname, paramarray);

    printf("finished output \n");
}

void getNGJitterMatrixEpochs(pulsar* pulse, int& NumEpochs)
{

    // count ECORR values
    if (pulse->nTNECORR > 0) {
        for (int i = 0; i < pulse->nTNECORR; i++) {
            printf("\nIncluding ECORR value for backend %s: %g mus", pulse->TNECORRFlagVal[i],
                   pulse->TNECORRVal[i]);
        }
    }

    // find number of epochs (default dt= 10 s)
    int* Processed = new int[pulse->nobs];

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    double dt = 10.0 / SECDAY;
    double satmin;
    double satmax;
    int nepoch = 0;
    int in = 0;
    int allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                in++;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }

    if (nepoch > 0) {
        printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
    }

    NumEpochs = nepoch;
}

void getNGSJitterMatrixEpochs(pulsar* pulse, int& NumEpochs)
{

    // count SECORR values
    if (pulse->nTNSECORR > 0) {
        for (int i = 0; i < pulse->nTNSECORR; i++) {
            printf("\nIncluding SECORR value for backend %s: %g mus", pulse->TNSECORRFlagVal[i],
                   pulse->TNSECORRVal[i]);
        }
    }

    // find number of epochs (default dt= 10 s)
    int* Processed = new int[pulse->nobs];

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNSECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    double dt = 10.0 / SECDAY;
    double satmin;
    double satmax;
    int nepoch = 0;
    int in = 0;
    int allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNSECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                in++;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }

    if (nepoch > 0) {
        printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
    }

    NumEpochs = nepoch;
}

void getNGJitterMatrix(pulsar* pulse, double** JitterMatrix, int& NumEpochs)
{

    // count ECORR values
    if (pulse->nTNECORR > 0) {
        for (int i = 0; i < pulse->nTNECORR; i++) {
            printf("\nIncluding ECORR value for backend %s: %g mus", pulse->TNECORRFlagVal[i],
                   pulse->TNECORRVal[i]);
        }
    }

    // find number of epochs (default dt= 10 s)
    int* Processed = new int[pulse->nobs];

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    double dt = 10.0 / SECDAY;
    double satmin;
    double satmax;
    int nepoch = 0;
    int in = 0;
    int allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                in++;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }

    if (nepoch > 0) {
        printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
    }

    NumEpochs = nepoch;

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    nepoch = 0;
    in = 0;
    allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                JitterMatrix[i][nepoch] = 1;
                                in++;
                            } else {
                                JitterMatrix[i][nepoch] = 0;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }
}

void getNGSJitterMatrix(pulsar* pulse, double** JitterMatrix, int& NumEpochs)
{

    // count SECORR values
    if (pulse->nTNSECORR > 0) {
        for (int i = 0; i < pulse->nTNSECORR; i++) {
            printf("\nIncluding SECORR value for backend %s: %g mus", pulse->TNSECORRFlagVal[i],
                   pulse->TNSECORRVal[i]);
        }
    }

    // find number of epochs (default dt= 10 s)
    int* Processed = new int[pulse->nobs];

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNSECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    double dt = 10.0 / SECDAY;
    double satmin;
    double satmax;
    int nepoch = 0;
    int in = 0;
    int allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNSECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                in++;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }

    if (nepoch > 0) {
        printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
    }

    NumEpochs = nepoch;

    // initialize processed flags
    for (int i = 0; i < pulse->nobs; i++) {
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i = 0; i < pulse->nobs; i++) {
        for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
            for (int k = 0; k < pulse->nTNSECORR; k++) {
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                    if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    nepoch = 0;
    in = 0;
    allProcessed = 0;
    while (!allProcessed) {
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i = 0; i < pulse->nobs; i++) {
            for (int j = 0; j < pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < pulse->nTNSECORR; k++) {
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNSECORRFlagID[k]) == 0) {
                        if (strcmp(pulse->obsn[i].flagVal[j], pulse->TNSECORRFlagVal[k]) == 0) {
                            if ((double)pulse->obsn[i].bat > satmin &&
                                (double)pulse->obsn[i].bat < satmax) {
                                Processed[i] = 1;
                                JitterMatrix[i][nepoch] = 1. / sqrt(pulse->obsn[i].tobs / 3600.);
                                in++;
                            } else {
                                JitterMatrix[i][nepoch] = 0;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0) {
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i = 0; i < pulse->nobs; i++) {
            if (Processed[i] == 0) {
                allProcessed = 0;
                break;
            }
        }
    }
}

void getCustomDMatrix(pulsar* pulse, int* MarginList, int** TempoFitNums, int* TempoJumpNums,
                      double** Dpriors, int incDM, int TimetoFit, int JumptoFit)
{

    double pdParamDeriv[MAX_PARAMS], dMultiplication;

    // Unset all fit flags for parameters we arn't marginalising over so they arn't in the design
    // Matrix

    int pcount = 1;
    int numToMargin = 1;

    Dpriors[0][0] = 0;
    Dpriors[0][1] = 0;

    for (int p = 1; p < TimetoFit; p++) {
        if (MarginList[pcount] != 1) {
            pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
        } else if (MarginList[pcount] == 1) {
            pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
            Dpriors[pcount][0] = 0;
            Dpriors[pcount][1] = 0;
            numToMargin++;
        }
        pcount++;
    }

    for (int i = 0; i < JumptoFit; i++) {
        if (MarginList[pcount] != 1) {
            pulse[0].fitJump[TempoJumpNums[i]] = 0;
        } else if (MarginList[pcount] == 1) {
            pulse[0].fitJump[TempoJumpNums[i]] = 1;
            Dpriors[pcount][0] = 0;
            Dpriors[pcount][1] = 0;
            numToMargin++;
        }
        pcount++;
    }

    //		for(int i=0; i < pulse->nobs; i++) {
    //			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv,
    // numToMargin, pulse, i,0); 			for(int j=0; j<numToMargin; j++) {
    // TNDM[i][j]=pdParamDeriv[j];
    //					//printf("Dmatrix: %i %i %22.20e \n", i,j,pdParamDeriv[j]);
    //			}
    //		}

    // Now set fit flags back to how they were

    for (int p = 1; p < TimetoFit; p++) {
        pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
    }

    for (int i = 0; i < JumptoFit; i++) {
        pulse[0].fitJump[TempoJumpNums[i]] = 1;
    }
}

void getEigenDVectorLike(void* context, Eigen::MatrixXd& TNDM, int Nobs, int TimeToMargin,
                         int TotalSize)
{

    int pcount = 0;
    int imargin = 0;

    pulsar* psr = ((MNStruct*)context)->pulse;
    FitInfo* fitinfo = &(psr->fitinfo);

    // we have to loop over parameters first then jumps
    // being careful to keep track of where we are.
    // This is because temponest keeps jumps at the end
    // but otherwise the parameters are in the same order.
    //
    // In tempo2 the jumps come before most parameters so we
    // skip over the part with the jumps without incrementing
    // the temponest index. Later we start from the total number
    // of parameters excluding jumps and only loop over the jumps

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        // this is something we want to marginalise over
        param_label p = fitinfo->paramIndex[iparam];
        // skip over parameters that are part of the TN model
        if (p == param_red_sin)
            continue;
        if (p == param_red_cos)
            continue;
        if (p == param_jitter)
            continue;
        if (p == param_red_dm_sin)
            continue;
        if (p == param_red_dm_cos)
            continue;

        // skip over jumps here. Note that jumps are at the start, so we skip without incrementing
        // pcount.
        if (p == param_JUMP)
            continue;

        if (((MNStruct*)context)->LDpriors[pcount][2] == 1) {

            const int k = fitinfo->paramCounters[iparam];
            for (int iobs = 0; iobs < psr->nobs; ++iobs) {
                const double x = psr->obsn[iobs].bat - psr->param[param_pepoch].val[0];

                TNDM(iobs, imargin) = fitinfo->paramDerivs[iparam](psr, 0, x, iobs, p, k);
            }
            ++imargin;
        }
        ++pcount;
    }

    // temponest has jumps at the end, so do NOT reset pcount!
    // pcount should be positioned at the first jump

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        param_label p = fitinfo->paramIndex[iparam];
        if (p == param_red_sin)
            continue;
        if (p == param_red_cos)
            continue;
        if (p == param_jitter)
            continue;
        if (p == param_red_dm_sin)
            continue;
        if (p == param_red_dm_cos)
            continue;

        if (p == param_JUMP) {
            if (((MNStruct*)context)->LDpriors[pcount][2] == 1) {

                const int k = fitinfo->paramCounters[iparam];
                for (int iobs = 0; iobs < psr->nobs; ++iobs) {
                    const double x = psr->obsn[iobs].bat - psr->param[param_pepoch].val[0];

                    TNDM(iobs, imargin) = fitinfo->paramDerivs[iparam](psr, 0, x, iobs, p, k);
                }

                ++imargin;
            }
            ++pcount;
        }
    }
}

void getCustomDVectorLike(void* context, double* TNDM, int Nobs, int TimeToMargin, int TotalSize)
{

    int pcount = 0;
    int imargin = 0;

    pulsar* psr = ((MNStruct*)context)->pulse;
    FitInfo* fitinfo = &(psr->fitinfo);

    // we have to loop over parameters first then jumps
    // being careful to keep track of where we are.
    // This is because temponest keeps jumps at the end
    // but otherwise the parameters are in the same order.
    //
    // In tempo2 the jumps come before most parameters so we
    // skip over the part with the jumps without incrementing
    // the temponest index. Later we start from the total number
    // of parameters excluding jumps and only loop over the jumps

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        // this is something we want to marginalise over
        param_label p = fitinfo->paramIndex[iparam];
        // skip over parameters that are part of the TN model
        if (p == param_red_sin)
            continue;
        if (p == param_red_cos)
            continue;
        if (p == param_jitter)
            continue;
        if (p == param_red_dm_sin)
            continue;
        if (p == param_red_dm_cos)
            continue;

        // skip over jumps here. Note that jumps are at the start, so we skip without incrementing
        // pcount.
        if (p == param_JUMP)
            continue;

        if (((MNStruct*)context)->LDpriors[pcount][2] == 1) {

            const int k = fitinfo->paramCounters[iparam];
            for (int iobs = 0; iobs < psr->nobs; ++iobs) {
                const double x = psr->obsn[iobs].bat - psr->param[param_pepoch].val[0];

                TNDM[iobs + imargin * ((MNStruct*)context)->pulse->nobs] =
                    fitinfo->paramDerivs[iparam](psr, 0, x, iobs, p, k);
                /// mjk - comment this out to make the code a lot faster
                // printf("TNDM: %i %i %i %g \n", iobs, imargin, iobs+imargin*((MNStruct
                // *)context)->pulse->nobs, fitinfo->paramDerivs[iparam](psr,0,x,iobs,p,k));
            }
            ++imargin;
        }
        ++pcount;
    }

    // temponest has jumps at the end, so do NOT reset pcount!
    // pcount should be positioned at the first jump

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        param_label p = fitinfo->paramIndex[iparam];
        if (p == param_red_sin)
            continue;
        if (p == param_red_cos)
            continue;
        if (p == param_jitter)
            continue;
        if (p == param_red_dm_sin)
            continue;
        if (p == param_red_dm_cos)
            continue;

        if (p == param_JUMP) {
            if (((MNStruct*)context)->LDpriors[pcount][2] == 1) {

                const int k = fitinfo->paramCounters[iparam];
                for (int iobs = 0; iobs < psr->nobs; ++iobs) {
                    const double x = psr->obsn[iobs].bat - psr->param[param_pepoch].val[0];

                    TNDM[iobs + imargin * ((MNStruct*)context)->pulse->nobs] =
                        fitinfo->paramDerivs[iparam](psr, 0, x, iobs, p, k);
                }

                ++imargin;
            }
            ++pcount;
        }
    }
}

void StoreTMatrix(double* TotalMatrix, void* context)
{

    int totalsize = ((MNStruct*)context)->totalsize;

    for (int i = 0; i < ((MNStruct*)context)->pulse->nobs * totalsize; i++) {
        TotalMatrix[i] = 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Form the Design Matrix////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    int TimetoMargin = ((MNStruct*)context)->TimetoMargin;
    if (TimetoMargin > 0) {
        Eigen::MatrixXd DMatrix(((MNStruct*)context)->pulse->nobs, TimetoMargin);
        getEigenDVectorLike(context, DMatrix, ((MNStruct*)context)->pulse->nobs, TimetoMargin,
                            totalsize);

        // Perform SVD
        Eigen::BDCSVD<Eigen::MatrixXd> svd(DMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::MatrixXd U = svd.matrixU();
        for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
            for (int j = 0; j < TimetoMargin; j++) {
                TotalMatrix[i + j * ((MNStruct*)context)->pulse->nobs] = U(i, j);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    double maxtspan = ((MNStruct*)context)->Tspan;

    int FitRedCoeff = 2 * (((MNStruct*)context)->numFitRedCoeff);
    int FitDMCoeff = 2 * (((MNStruct*)context)->numFitDMCoeff);
    int FitScatCoeff = 2 * (((MNStruct*)context)->numFitScatCoeff);
    int FitBandCoeff = 2 * (((MNStruct*)context)->numFitBandNoiseCoeff);
    int FitGroupNoiseCoeff = 2 * ((MNStruct*)context)->numFitGroupNoiseCoeff;

    int totCoeff = ((MNStruct*)context)->totCoeff;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double* freqs = new double[totCoeff];
    double* DMVec = new double[((MNStruct*)context)->pulse->nobs];

    double DMKappa = 2.410 * std::pow(10.0, -16);
    int startpos = 0;

    if (((MNStruct*)context)->incRED > 0 || ((MNStruct*)context)->incGWB == 1) {
        for (int i = 0; i < FitRedCoeff / 2 - ((MNStruct*)context)->incFloatRed; i++) {

            freqs[startpos + i] = (double)((MNStruct*)context)->sampleFreq[i] / maxtspan;
            freqs[startpos + i + FitRedCoeff / 2] = freqs[startpos + i];
        }

        for (int i = 0; i < FitRedCoeff / 2; i++) {
            for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
                double time = (double)((MNStruct*)context)->pulse->obsn[k].bat;
                TotalMatrix[k + (i + TimetoMargin + startpos) * ((MNStruct*)context)->pulse->nobs] =
                    cos(2 * M_PI * freqs[i] * time);
                TotalMatrix[k + (i + FitRedCoeff / 2 + TimetoMargin + startpos) *
                                    ((MNStruct*)context)->pulse->nobs] =
                    sin(2 * M_PI * freqs[i] * time);
            }
        }

        startpos += FitRedCoeff;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)context)->incDM > 0) {

        for (int o = 0; o < ((MNStruct*)context)->pulse->nobs; o++) {
            DMVec[o] =
                1.0 / (DMKappa * std::pow((double)((MNStruct*)context)->pulse->obsn[o].freqSSB, 2));
        }

        for (int i = 0; i < FitDMCoeff / 2; i++) {

            freqs[startpos + i] =
                ((MNStruct*)context)
                    ->sampleFreq[startpos / 2 - ((MNStruct*)context)->incFloatRed + i] /
                maxtspan;
            freqs[startpos + i + FitDMCoeff / 2] = freqs[startpos + i];

            for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
                double time = (double)((MNStruct*)context)->pulse->obsn[k].bat;

                TotalMatrix[k + (i + TimetoMargin + startpos) * ((MNStruct*)context)->pulse->nobs] =
                    cos(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
                TotalMatrix[k + (i + FitDMCoeff / 2 + TimetoMargin + startpos) *
                                    ((MNStruct*)context)->pulse->nobs] =
                    sin(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
            }
        }

        startpos += FitDMCoeff;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Scattering variations ////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)context)->incScat > 0) {
        double* ScatVec = new double[((MNStruct*)context)->pulse->nobs];
        for (int o = 0; o < ((MNStruct*)context)->pulse->nobs; o++)
            ScatVec[o] = 1. / (pow((double)((MNStruct*)context)->pulse->obsn[o].freqSSB / 1400e6,
                                   4));  // Referenced to 1.4 Ghz to be consistent with Enterprise;
                                         // Invert to be consistent with Tempo2 - 20220810 GD/AP

        for (int i = 0; i < FitScatCoeff / 2; i++) {
            freqs[startpos + i] =
                ((MNStruct*)context)
                    ->sampleFreq[startpos / 2 - ((MNStruct*)context)->incFloatRed + i] /
                maxtspan;
            freqs[startpos + i + FitScatCoeff / 2] = freqs[startpos + i];

            for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
                double time = (double)((MNStruct*)context)->pulse->obsn[k].bat;
                TotalMatrix[k + (i + TimetoMargin + startpos) * ((MNStruct*)context)->pulse->nobs] =
                    cos(2 * M_PI * freqs[startpos + i] * time) * ScatVec[k];
                TotalMatrix[k + (i + FitScatCoeff / 2 + TimetoMargin + startpos) *
                                    ((MNStruct*)context)->pulse->nobs] =
                    sin(2 * M_PI * freqs[startpos + i] * time) * ScatVec[k];
                // printf("tot Mat = %lg\n", TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct
                // *)context)->pulse->nobs]);
            }
        }
        startpos += FitScatCoeff;
        delete[] ScatVec;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Band DM/////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)context)->incBandNoise > 0) {

        if (((MNStruct*)context)->incDM == 0) {
            for (int o = 0; o < ((MNStruct*)context)->pulse->nobs; o++) {
                DMVec[o] =
                    1.0 /
                    (DMKappa * std::pow((double)((MNStruct*)context)->pulse->obsn[o].freqSSB, 2));
            }
        }

        for (int b = 0; b < ((MNStruct*)context)->incBandNoise; b++) {

            double startfreq = ((MNStruct*)context)->FitForBand[b][0];
            double stopfreq = ((MNStruct*)context)->FitForBand[b][1];
            double BandScale = ((MNStruct*)context)->FitForBand[b][2];
            int BandPriorType = ((MNStruct*)context)->FitForBand[b][3];

            double Tspan = maxtspan;
            double f1yr = 1.0 / 3.16e7;

            for (int i = 0; i < FitBandCoeff / 2; i++) {

                freqs[startpos + i] = ((double)(i + 1)) / maxtspan;
                freqs[startpos + i + FitBandCoeff / 2] = freqs[startpos + i];
            }

            for (int i = 0; i < FitBandCoeff / 2; i++) {
                for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
                    if (((MNStruct*)context)->pulse->obsn[k].freq > startfreq &&
                        ((MNStruct*)context)->pulse->obsn[k].freq < stopfreq) {
                        double time = (double)((MNStruct*)context)->pulse->obsn[k].bat;
                        TotalMatrix[k + (i + TimetoMargin + startpos) *
                                            ((MNStruct*)context)->pulse->nobs] =
                            cos(2 * M_PI * freqs[startpos + i] * time);
                        TotalMatrix[k + (i + TimetoMargin + startpos + FitBandCoeff / 2) *
                                            ((MNStruct*)context)->pulse->nobs] =
                            sin(2 * M_PI * freqs[startpos + i] * time);
                    } else {
                        TotalMatrix[k + (i + TimetoMargin + startpos) *
                                            ((MNStruct*)context)->pulse->nobs] = 0;
                        TotalMatrix[k + (i + TimetoMargin + startpos + FitBandCoeff / 2) *
                                            ((MNStruct*)context)->pulse->nobs] = 0;
                    }
                }
            }

            startpos += FitBandCoeff;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Add Group Noise/////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)context)->incGroupNoise > 0) {

        for (int g = 0; g < ((MNStruct*)context)->incGroupNoise; g++) {

            startpos = startpos + FitGroupNoiseCoeff;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Add ECORR Coeffs////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (((MNStruct*)context)->incNGJitter > 0) {

        // printf("Calling\n");
        double** NGJitterMatrix;
        int NumNGEpochs = ((MNStruct*)context)->numNGJitterEpochs;
        // getNGJitterMatrixEpochs(psr, NumNGEpochs);

        NGJitterMatrix = new double*[((MNStruct*)context)->pulse->nobs];
        for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
            NGJitterMatrix[i] = new double[NumNGEpochs];
            for (int j = 0; j < NumNGEpochs; j++) {
                NGJitterMatrix[i][j] = 0;
            }
        }
        int* NGJitterSysFlags = new int[((MNStruct*)context)->pulse->nobs];

        for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
            for (int j = 0; j < ((MNStruct*)context)->pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < ((MNStruct*)context)->pulse->nTNECORR; k++) {
                    if (strcmp(((MNStruct*)context)->pulse->obsn[i].flagID[j],
                               ((MNStruct*)context)->pulse->TNECORRFlagID[k]) == 0) {
                        if (strcmp(((MNStruct*)context)->pulse->obsn[i].flagVal[j],
                                   ((MNStruct*)context)->pulse->TNECORRFlagVal[k]) == 0) {
                            NGJitterSysFlags[i] = k;
                            // printf("NGFLag? %i %i %s \n", i, k, psr->TNECORRFlagVal[k]);
                        }
                    }
                }
            }
        }

        getNGJitterMatrix(((MNStruct*)context)->pulse, NGJitterMatrix, NumNGEpochs);

        int* NGJitterEpochFlags = new int[NumNGEpochs];

        for (int i = 0; i < NumNGEpochs; i++) {
            for (int j = 0; j < ((MNStruct*)context)->pulse->nobs; j++) {
                if (NGJitterMatrix[j][i] != 0) {
                    NGJitterEpochFlags[i] = NGJitterSysFlags[j];
                }
            }
        }

        ((MNStruct*)context)->NGJitterSysFlags = NGJitterEpochFlags;

        for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
            for (int i = 0; i < ((MNStruct*)context)->numNGJitterEpochs; i++) {
                TotalMatrix[k + (i + TimetoMargin + startpos) * ((MNStruct*)context)->pulse->nobs] =
                    NGJitterMatrix[k][i];  /// sqrt((((MNStruct *)context)->TobsInfo[k]/3600.0));
            }
        }

        delete[] NGJitterSysFlags;
        for (int j = 0; j < ((MNStruct*)context)->pulse->nobs; j++) {
            delete[] NGJitterMatrix[j];
        }
        delete[] NGJitterMatrix;
    }

    //// ADDING SECORR

    if (((MNStruct*)context)->incNGSJitter > 0) {

        // printf("Calling\n");
        double** NGSJitterMatrix;
        int NumNGSEpochs = ((MNStruct*)context)->numNGSJitterEpochs;
        // getNGJitterMatrixEpochs(psr, NumNGEpochs);

        NGSJitterMatrix = new double*[((MNStruct*)context)->pulse->nobs];
        for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
            NGSJitterMatrix[i] = new double[NumNGSEpochs];
            for (int j = 0; j < NumNGSEpochs; j++) {
                NGSJitterMatrix[i][j] = 0;
            }
        }
        int* NGSJitterSysFlags = new int[((MNStruct*)context)->pulse->nobs];

        for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
            for (int j = 0; j < ((MNStruct*)context)->pulse->obsn[i].nFlags; j++) {
                for (int k = 0; k < ((MNStruct*)context)->pulse->nTNSECORR; k++) {
                    if (strcmp(((MNStruct*)context)->pulse->obsn[i].flagID[j],
                               ((MNStruct*)context)->pulse->TNSECORRFlagID[k]) == 0) {
                        if (strcmp(((MNStruct*)context)->pulse->obsn[i].flagVal[j],
                                   ((MNStruct*)context)->pulse->TNSECORRFlagVal[k]) == 0) {
                            NGSJitterSysFlags[i] = k;
                            // printf("NGFLag? %i %i %s \n", i, k, psr->TNECORRFlagVal[k]);
                        }
                    }
                }
            }
        }

        getNGSJitterMatrix(((MNStruct*)context)->pulse, NGSJitterMatrix, NumNGSEpochs);

        int* NGSJitterEpochFlags = new int[NumNGSEpochs];

        for (int i = 0; i < NumNGSEpochs; i++) {
            for (int j = 0; j < ((MNStruct*)context)->pulse->nobs; j++) {
                if (NGSJitterMatrix[j][i] != 0) {
                    NGSJitterEpochFlags[i] = NGSJitterSysFlags[j];
                }
            }
        }

        ((MNStruct*)context)->NGSJitterSysFlags = NGSJitterEpochFlags;

        for (int k = 0; k < ((MNStruct*)context)->pulse->nobs; k++) {
            for (int i = 0; i < ((MNStruct*)context)->numNGSJitterEpochs; i++) {
                TotalMatrix[k + (i + TimetoMargin + startpos) * ((MNStruct*)context)->pulse->nobs] =
                    NGSJitterMatrix[k][i];  /// sqrt((((MNStruct *)context)->TobsInfo[k]/3600.0));
            }
        }

        delete[] NGSJitterSysFlags;
        for (int j = 0; j < ((MNStruct*)context)->pulse->nobs; j++) {
            delete[] NGSJitterMatrix[j];
        }
        delete[] NGSJitterMatrix;
    }

    delete[] DMVec;
    delete[] freqs;
}

void getArraySizeInfo(void* context)
{

    int TimetoMargin = 0;
    for (int i = 0; i < ((MNStruct*)context)->numFitTiming + ((MNStruct*)context)->numFitJumps;
         i++) {
        if (((MNStruct*)context)->LDpriors[i][2] == 1)
            TimetoMargin++;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double start, end;
    int go = 0;
    for (int i = 0; i < ((MNStruct*)context)->pulse->nobs; i++) {
        if (((MNStruct*)context)->pulse->obsn[i].deleted == 0) {
            if (go == 0) {
                go = 1;
                start = (double)((MNStruct*)context)->pulse->obsn[i].bat;
                end = start;
            } else {
                if (start > (double)((MNStruct*)context)->pulse->obsn[i].bat)
                    start = (double)((MNStruct*)context)->pulse->obsn[i].bat;
                if (end < (double)((MNStruct*)context)->pulse->obsn[i].bat)
                    end = (double)((MNStruct*)context)->pulse->obsn[i].bat;
            }
        }
    }

    double maxtspan = 1 * (end - start);

    int FitRedCoeff = 2 * (((MNStruct*)context)->numFitRedCoeff);
    int FitDMCoeff = 2 * (((MNStruct*)context)->numFitDMCoeff);
    int FitScatCoeff = 2 * (((MNStruct*)context)->numFitScatCoeff);
    int FitBandNoiseCoeff = 2 * (((MNStruct*)context)->numFitBandNoiseCoeff);
    int FitGroupNoiseCoeff = 2 * ((MNStruct*)context)->numFitGroupNoiseCoeff;

    int NumNGEpochs = 0;
    if (((MNStruct*)context)->incNGJitter > 0) {
        getNGJitterMatrixEpochs(((MNStruct*)context)->pulse, NumNGEpochs);
    }
    ((MNStruct*)context)->numNGJitterEpochs = NumNGEpochs;

    int NumNGSEpochs = 0;
    if (((MNStruct*)context)->incNGSJitter > 0) {
        getNGSJitterMatrixEpochs(((MNStruct*)context)->pulse, NumNGSEpochs);
    }
    ((MNStruct*)context)->numNGSJitterEpochs = NumNGSEpochs;

    int totCoeff = 0;
    if (((MNStruct*)context)->incRED != 0 || ((MNStruct*)context)->incGWB == 1)
        totCoeff += FitRedCoeff;
    if (((MNStruct*)context)->incDM != 0)
        totCoeff += FitDMCoeff;
    if (((MNStruct*)context)->incScat != 0)
        totCoeff += FitScatCoeff;
    if (((MNStruct*)context)->incBandNoise > 0)
        totCoeff += ((MNStruct*)context)->incBandNoise * FitBandNoiseCoeff;
    if (((MNStruct*)context)->incNGJitter > 0)
        totCoeff += ((MNStruct*)context)->numNGJitterEpochs;
    if (((MNStruct*)context)->incNGSJitter > 0)
        totCoeff += ((MNStruct*)context)->numNGSJitterEpochs;
    if (((MNStruct*)context)->incGroupNoise > 0)
        totCoeff += ((MNStruct*)context)->incGroupNoise * FitGroupNoiseCoeff;

    int MarginRedShapelets = ((MNStruct*)context)->MarginRedShapeCoeff;
    int totalredshapecoeff = 0;

    if (((MNStruct*)context)->incRedShapeEvent != 0) {
        if (MarginRedShapelets == 1) {
            totalredshapecoeff =
                ((MNStruct*)context)->numRedShapeCoeff * ((MNStruct*)context)->incRedShapeEvent;
        }
    }

    int totalsize = TimetoMargin + totCoeff + totalredshapecoeff;

    ((MNStruct*)context)->Tspan = maxtspan;
    ((MNStruct*)context)->TimetoMargin = TimetoMargin;
    ((MNStruct*)context)->totCoeff = totCoeff;
    ((MNStruct*)context)->totRedShapeCoeff = totalredshapecoeff;
    ((MNStruct*)context)->totalsize = totalsize;

    //	printf("TimetoMargin %i, totCoeff %i, totalredshapecoeff %i, totalsize %i \n", TimetoMargin,
    // totCoeff, totalredshapecoeff, totalsize);
}
