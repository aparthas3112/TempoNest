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
#include "eigen_config.h"
#include "tempo2.h"
#include "types/model.h"

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
}

void readsummary(pulsar* psr, std::string longname, int ndim, void* context, int ndims)
{

    std::vector<double> paramlist(2 * ndims);

    double** paramarray = new double*[ndims];
    for (int p = 0; p < ndims; p++) {
        paramarray[p] = new double[4];
    }

    readtxtoutput(longname, ndim, paramarray);
    readphyslive(longname, ndim, paramarray, 0);

    formBatsAll(globals::pulsar, 1);  // Form Barycentric arrival times
    // printf("formed bats \n");
    formResiduals(globals::pulsar, 1, 1);  // Form residuals

    double Evidence = 0;
    TNtextOutput(globals::pulsar, 1, 0, context, ndims, paramlist, Evidence, longname, paramarray);

    printf("finished output \n");
}

void getEigenDVectorLike(Eigen::MatrixXd& design_matrix)
{

    timing_model_t* timing_model = model::timing_model->as<timing_model_t>();

    int imargin = 0;

    FitInfo* fitinfo = &(globals::pulsar->fitinfo);

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        // check if this is something we want to marginalise over
        if (!timing_model->marginalised[iparam]) {
            continue;
        }

        param_label p = fitinfo->paramIndex[iparam];

        const int k = fitinfo->paramCounters[iparam];
        for (int iobs = 0; iobs < globals::pulsar->nobs; ++iobs) {
            const double x =
                globals::pulsar->obsn[iobs].bat - globals::pulsar->param[param_pepoch].val[0];

            design_matrix(iobs, imargin) =
                fitinfo->paramDerivs[iparam](globals::pulsar, 0, x, iobs, p, k);
        }
        ++imargin;
    }
}

void StoreTMatrix()
{

    int totalsize = model::total_size;

    model::total_matrix = Eigen::MatrixXd::Zero(globals::pulsar->nobs, totalsize);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Form the Design Matrix////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Forming Design Matrix " << model::design_size << std::endl;
    int TimetoMargin = model::design_size;
    if (TimetoMargin > 0) {
        model::design_matrix = Eigen::MatrixXd::Zero(globals::pulsar->nobs, TimetoMargin);
        getEigenDVectorLike(model::design_matrix);

        // Perform SVD
        Eigen::BDCSVD<Eigen::MatrixXd> svd(model::design_matrix,
                                           Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::MatrixXd U = svd.matrixU();
        for (int i = 0; i < globals::pulsar->nobs; i++) {
            for (int j = 0; j < TimetoMargin; j++) {
                model::total_matrix(i, j) = U(i, j);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    double maxtspan = model::max_tspan;
    int totCoeff = model::noise_size;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double* freqs = new double[totCoeff];
    double* DMVec = new double[globals::pulsar->nobs];

    double DMKappa = 2.410 * std::pow(10.0, -16);
    int startpos = 0;

    if (model::pl_red_noise.has_value()) {
        pl_red_noise_element* pl = model::pl_red_noise.value()->as<pl_red_noise_element>();

        for (int i = 0; i < pl->num_freqs; i++) {

            freqs[startpos + i] = pl->frequencies[i] / maxtspan;
            freqs[startpos + i + pl->num_freqs] = freqs[startpos + i];
        }

        for (int i = 0; i < pl->num_freqs; i++) {
            for (int k = 0; k < globals::pulsar->nobs; k++) {
                double time = (double)globals::pulsar->obsn[k].bat;
                model::total_matrix(k, i + TimetoMargin + startpos) =
                    cos(2 * M_PI * freqs[i] * time);
                model::total_matrix(k, i + pl->num_freqs + TimetoMargin + startpos) =
                    sin(2 * M_PI * freqs[i] * time);
            }
        }

        startpos += 2 * pl->num_freqs;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (model::pl_dm_noise.has_value()) {

        pl_dm_noise_element* pl = model::pl_dm_noise.value()->as<pl_dm_noise_element>();

        for (int o = 0; o < globals::pulsar->nobs; o++) {
            DMVec[o] = 1.0 / (DMKappa * std::pow((double)globals::pulsar->obsn[o].freqSSB, 2));
        }

        for (int i = 0; i < pl->num_freqs; i++) {

            freqs[startpos + i] = pl->frequencies[i] / maxtspan;
            freqs[startpos + i + pl->num_freqs] = freqs[startpos + i];

            for (int k = 0; k < globals::pulsar->nobs; k++) {
                double time = (double)globals::pulsar->obsn[k].bat;

                model::total_matrix(k, i + TimetoMargin + startpos) =
                    cos(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];

                model::total_matrix(k, i + pl->num_freqs + TimetoMargin + startpos) =
                    sin(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
            }
        }

        startpos += 2 * pl->num_freqs;
    }

    delete[] DMVec;
    delete[] freqs;
}

void getArraySizeInfo()
{

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    timing_model_t* timing_model = model::timing_model->as<timing_model_t>();
    model::design_size = timing_model->design_size;

    double start, end;
    int go = 0;
    for (int i = 0; i < globals::pulsar->nobs; i++) {
        if (globals::pulsar->obsn[i].deleted == 0) {
            if (go == 0) {
                go = 1;
                start = (double)globals::pulsar->obsn[i].bat;
                end = start;
            } else {
                if (start > (double)globals::pulsar->obsn[i].bat)
                    start = (double)globals::pulsar->obsn[i].bat;
                if (end < (double)globals::pulsar->obsn[i].bat)
                    end = (double)globals::pulsar->obsn[i].bat;
            }
        }
    }

    double maxtspan = 1 * (end - start);

    int totCoeff = 0;
    if (model::pl_red_noise.has_value()) {
        pl_red_noise_element* pl = model::pl_red_noise.value()->as<pl_red_noise_element>();
        totCoeff += 2 * pl->num_freqs;
    }

    if (model::pl_dm_noise.has_value()) {
        pl_dm_noise_element* pl = model::pl_dm_noise.value()->as<pl_dm_noise_element>();
        totCoeff += 2 * pl->num_freqs;
    }

    int totalsize = model::design_size + totCoeff;

    model::max_tspan = maxtspan;
    model::noise_size = totCoeff;
    model::total_size = totalsize;
}
