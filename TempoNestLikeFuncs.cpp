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
 *    then please acknowledge it by citing Lentati L., Alexander P., Hobson M.
 * P. (2013) for TempoNest, Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369,
 * Issue 2, pp. 655-672 (bibtex: 2006MNRAS.369..655H) or Edwards, Hobbs &
 * Manchester (2006) MNRAS, VOl 372, Issue 4, pp. 1549-1574 (bibtex:
 * 2006MNRAS.372.1549E) when discussing the timing model and MultiNest Papers
 * here.
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
#include "dgemm.h"
#include "dgemv.h"
#include "dgesvd.h"
#include "dpotrf.h"
#include "dpotri.h"
#include "dpotrs.h"
#include "qrdecomp.h"
#include "tempo2.h"

#ifdef HAVE_MLAPACK
#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>
#endif
#include <iostream>

using namespace std;

void *globalcontext;
// void SmallNelderMeadOptimum(int nParameters, double *ML);
// void  WriteProfileFreqEvo(std::string longname, int &ndim, int
// profiledimstart);

void assigncontext(void *context) { globalcontext = context; }

/*
int Wrap(int kX, int const kLowerBound, int const kUpperBound)
{
    int range_size = kUpperBound - kLowerBound + 1;

    if (kX < kLowerBound)
        kX += range_size * ((kLowerBound - kX) / range_size + 1);

    return kLowerBound + (kX - kLowerBound) % range_size;
}

*/
void TNothpl(int n, double x, double *pl) {
  double a = 2.0;
  double b = 0.0;
  double c = 1.0;
  double y0 = 1.0;
  double y1 = 2.0 * x;
  pl[0] = 1.0;
  pl[1] = 2.0 * x;

  //	printf("I AM IN OTHPL %i \n", n);
  for (int k = 2; k < n; k++) {
    double c = 2.0 * (k - 1.0);
    //		printf("%i %g\n", k, sqrt(double(k*1.0)));
    y0 = y0 / std::sqrt(double(k * 1.0));
    y1 = y1 / std::sqrt(double(k * 1.0));
    double yn = (a * x + b) * y1 - c * y0;
    yn = yn;     /// sqrt(double(k));
    pl[k] = yn;  /// sqrt(double(k));
    y0 = y1;
    y1 = yn;
  }
}

void TNothplMC(int n, double x, double *pl, int cpos) {
  // printf("I AM IN OTHPL %i %i\n", n, cpos);
  double a = 2.0;
  double b = 0.0;
  double c = 1.0;
  double y0 = 1.0;
  double y1 = 2.0 * x;
  pl[0 + cpos] = 1.0;
  if (n > 1) {
    pl[1 + cpos] = 2.0 * x;
  }

  for (int k = 2; k < n; k++) {
    double c = 2.0 * (k - 1.0);
    //		printf("%i %g\n", k, sqrt(double(k*1.0)));
    y0 = y0 / std::sqrt(double(k * 1.0));
    y1 = y1 / std::sqrt(double(k * 1.0));
    double yn = (a * x + b) * y1 - c * y0;
    yn = yn;            /// sqrt(double(k));
    pl[k + cpos] = yn;  /// sqrt(double(k));
    y0 = y1;
    y1 = yn;
  }
}

double FastNewLRedMarginLogLike(double *Cube, int ndim, double *DerivedParams,
                                int npars, void *context) {
  printf("In fast like\n");
  int pcount = 0;

  int numfit = ((MNStruct *)globalcontext)->numFitTiming +
               ((MNStruct *)globalcontext)->numFitJumps;
  long double LDparams[numfit];
  double Fitparams[numfit];
  double *Resvec = new double[((MNStruct *)globalcontext)->pulse->nobs];
  int fitcount = 0;

  pcount = 0;

  for (int p = 0; p < ((MNStruct *)globalcontext)->numFitTiming +
                          ((MNStruct *)globalcontext)->numFitJumps;
       p++) {
    if (((MNStruct *)globalcontext)->Dpriors[p][1] !=
        ((MNStruct *)globalcontext)->Dpriors[p][0]) {
      double val = 0;
      if ((((MNStruct *)globalcontext)->LDpriors[p][3]) == 0) {
        val = Cube[fitcount];
      }
      if ((((MNStruct *)globalcontext)->LDpriors[p][3]) == 1) {
        val = pow(10.0, Cube[fitcount]);
      }
      LDparams[p] = val * (((MNStruct *)globalcontext)->LDpriors[p][1]) +
                    (((MNStruct *)globalcontext)->LDpriors[p][0]);
      fitcount++;
    } else if (((MNStruct *)globalcontext)->Dpriors[p][1] ==
               ((MNStruct *)globalcontext)->Dpriors[p][0]) {
      LDparams[p] = ((MNStruct *)globalcontext)->Dpriors[p][0] *
                        (((MNStruct *)globalcontext)->LDpriors[p][1]) +
                    (((MNStruct *)globalcontext)->LDpriors[p][0]);
    }
  }
  pcount = 0;
  double phase = (double)LDparams[0];
  pcount++;
  for (int p = 1; p < ((MNStruct *)globalcontext)->numFitTiming; p++) {
    ((MNStruct *)globalcontext)
        ->pulse->param[((MNStruct *)globalcontext)->TempoFitNums[p][0]]
        .val[((MNStruct *)globalcontext)->TempoFitNums[p][1]] =
        LDparams[pcount];
    pcount++;
  }
  for (int p = 0; p < ((MNStruct *)globalcontext)->numFitJumps; p++) {
    ((MNStruct *)globalcontext)
        ->pulse->jumpVal[((MNStruct *)globalcontext)->TempoJumpNums[p]] =
        LDparams[pcount];
    pcount++;
  }

  fastformBatsAll(((MNStruct *)globalcontext)->pulse,
                  ((MNStruct *)globalcontext)->numberpulsars);
  formResiduals(((MNStruct *)globalcontext)->pulse,
                ((MNStruct *)globalcontext)->numberpulsars, 1);

  for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
    Resvec[o] =
        (double)((MNStruct *)globalcontext)->pulse->obsn[o].residual + phase;
  }

  pcount = fitcount;

  if (((MNStruct *)globalcontext)->incStep > 0) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->incStep; i++) {
      double StepAmp = Cube[pcount];
      pcount++;
      double StepTime = Cube[pcount];
      pcount++;
      for (int o1 = 0; o1 < ((MNStruct *)globalcontext)->pulse->nobs; o1++) {
        if (((MNStruct *)globalcontext)->pulse->obsn[o1].bat > StepTime) {
          Resvec[o1] += StepAmp;
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Subtract
  /// GLitches/////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < ((MNStruct *)globalcontext)->incGlitch; i++) {
    double GlitchMJD = Cube[pcount];
    pcount++;

    double *GlitchAmps = new double[2];
    if (((MNStruct *)globalcontext)->incGlitchTerms == 1 ||
        ((MNStruct *)globalcontext)->incGlitchTerms == 3) {
      GlitchAmps[0] = Cube[pcount];
      pcount++;
    } else if (((MNStruct *)globalcontext)->incGlitchTerms == 2) {
      GlitchAmps[0] = Cube[pcount];
      pcount++;
      GlitchAmps[1] = Cube[pcount];
      pcount++;
    }

    for (int o1 = 0; o1 < ((MNStruct *)globalcontext)->pulse->nobs; o1++) {
      if (((MNStruct *)globalcontext)->pulse->obsn[o1].bat > GlitchMJD) {
        if (((MNStruct *)globalcontext)->incGlitchTerms == 1) {
          long double arg = 0;
          arg =
              ((((MNStruct *)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
               ((MNStruct *)globalcontext)->pulse->param[param_f].val[0]) *
              86400.0;
          double darg = (double)arg;
          Resvec[o1] += GlitchAmps[0] * darg;

        } else if (((MNStruct *)globalcontext)->incGlitchTerms == 2) {
          for (int j = 0; j < ((MNStruct *)globalcontext)->incGlitchTerms;
               j++) {
            long double arg = 0;
            if (j == 0) {
              arg =
                  ((((MNStruct *)globalcontext)->pulse->obsn[o1].bat -
                    GlitchMJD) /
                   ((MNStruct *)globalcontext)->pulse->param[param_f].val[0]) *
                  86400.0;
            }
            if (j == 1) {
              arg = 0.5 *
                    pow((((MNStruct *)globalcontext)->pulse->obsn[o1].bat -
                         GlitchMJD) *
                            86400.0,
                        2) /
                    ((MNStruct *)globalcontext)->pulse->param[param_f].val[0];
            }
            double darg = (double)arg;
            Resvec[o1] += GlitchAmps[j] * darg;
          }
        } else if (((MNStruct *)globalcontext)->incGlitchTerms == 3) {
          long double arg = 0;
          arg = 0.5 *
                pow((((MNStruct *)globalcontext)->pulse->obsn[o1].bat -
                     GlitchMJD) *
                        86400.0,
                    2) /
                ((MNStruct *)globalcontext)->pulse->param[param_f].val[0];
          double darg = (double)arg;
          Resvec[o1] += GlitchAmps[0] * darg;
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Get White Noise
  /// vector///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double *Noise;
  Noise = new double[((MNStruct *)globalcontext)->pulse->nobs];

  for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
    Noise[o] = ((MNStruct *)globalcontext)->PreviousNoise[o];
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Get Time domain
  /// likelihood//////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double tdet = 0;
  double timelike = 0;

  for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
    timelike += Resvec[o] * Resvec[o] * Noise[o];
    tdet -= log(Noise[o]);
  }

  printf("Fast time like %g %g %i\n", timelike, tdet,
         ((MNStruct *)globalcontext)->totCoeff);

  //////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////Do
  /// Algebra/////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  int totCoeff = ((MNStruct *)globalcontext)->totCoeff;

  int TimetoMargin = 0;
  for (int i = 0; i < ((MNStruct *)globalcontext)->numFitTiming +
                          ((MNStruct *)globalcontext)->numFitJumps;
       i++) {
    if (((MNStruct *)globalcontext)->LDpriors[i][2] == 1) TimetoMargin++;
  }

  int totalsize = totCoeff + TimetoMargin;

  double *NTd = new double[totalsize];

  double **NT = new double *[((MNStruct *)globalcontext)->pulse->nobs];
  for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
    NT[i] = new double[totalsize];
    for (int j = 0; j < totalsize; j++) {
      NT[i][j] = ((MNStruct *)globalcontext)->PreviousNT[i][j];
    }
  }

  double **TNT = new double *[totalsize];
  for (int i = 0; i < totalsize; i++) {
    TNT[i] = new double[totalsize];
    for (int j = 0; j < totalsize; j++) {
      TNT[i][j] = ((MNStruct *)globalcontext)->PreviousTNT[i][j];
    }
  }

  dgemv(NT, Resvec, NTd, ((MNStruct *)globalcontext)->pulse->nobs, totalsize,
        'T');

  double freqlike = 0;
  double *WorkCoeff = new double[totalsize];
  for (int o1 = 0; o1 < totalsize; o1++) {
    WorkCoeff[o1] = NTd[o1];
  }

  int globalinfo = 0;
  int info = 0;
  double jointdet = ((MNStruct *)globalcontext)->PreviousJointDet;
  double freqdet = ((MNStruct *)globalcontext)->PreviousFreqDet;
  double uniformpriorterm = ((MNStruct *)globalcontext)->PreviousUniformPrior;

  info = 0;
  dpotrsInfo(TNT, WorkCoeff, totalsize, info);
  if (info != 0) globalinfo = 1;

  for (int j = 0; j < totalsize; j++) {
    freqlike += NTd[j] * WorkCoeff[j];
  }

  double lnew = -0.5 * (tdet + jointdet + freqdet + timelike - freqlike) +
                uniformpriorterm;
  std::cout << "lnew " << lnew << " " << tdet << " " << jointdet << " "
            << freqdet << " " << timelike << " " << freqlike << " "
            << uniformpriorterm << std::endl;

  if (isnan(lnew) || isinf(lnew) || globalinfo != 0) {
    lnew = -pow(10.0, 20);
  }

  /*
          if(((MNStruct *)globalcontext)->uselongdouble > 0 ){
                  fpu_fix_end(&oldcw);
          }
  */

  delete[] WorkCoeff;
  delete[] NTd;
  delete[] Noise;
  delete[] Resvec;

  for (int j = 0; j < ((MNStruct *)globalcontext)->pulse->nobs; j++) {
    delete[] NT[j];
  }
  delete[] NT;

  for (int j = 0; j < totalsize; j++) {
    delete[] TNT[j];
  }
  delete[] TNT;

  return lnew;
}

void LRedLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew,
                    void *context) {
  for (int p = 0; p < ndim; p++) {
    Cube[p] = (((MNStruct *)globalcontext)->PriorsArray[p + ndim] -
               ((MNStruct *)globalcontext)->PriorsArray[p]) *
                  Cube[p] +
              ((MNStruct *)globalcontext)->PriorsArray[p];
  }

  double *DerivedParams = new double[npars];

  // double result = NewLRedMarginLogLike(ndim, Cube, npars, DerivedParams,
  // context);
  double result =
      NewLRedMarginLogLike(Cube, ndim, DerivedParams, npars, context);

  delete[] DerivedParams;

  lnew = result;
}

// double  NewLRedMarginLogLike(int &ndim, double *Cube, int &npars, double
// *DerivedParams, void *context){
double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived,
                            void *context) {
  logtchk("Entering TempoNest likelihood");

  double uniformpriorterm = 0;
  clock_t startClock, endClock;

  double **EFAC;
  double *EQUAD;
  int pcount = 0;

  int numfit = ((MNStruct *)globalcontext)->numFitTiming +
               ((MNStruct *)globalcontext)->numFitJumps;
  int TimetoMargin = ((MNStruct *)globalcontext)->TimetoMargin;
  long double LDparams[numfit];
  for (int i = 0; i < numfit; i++) {
    LDparams[i] = 0;
  }
  double Fitparams[numfit];
  double *Resvec = new double[((MNStruct *)globalcontext)->pulse->nobs];
  int fitcount = 0;

  pcount = 0;
  ///////////////hacky glitch thing for Vela, comment this out for *anything
  /// else*/////////////////////////
  // for(int p=1;p<9; p++){
  // Cube[p] = Cube[0];
  //	printf("Cube 0 %g \n", Cube[0]);
  //}
  // for(int p=10;p<17; p++){
  //      Cube[p] = Cube[9];
  //	printf("Cube 9 %g\n", Cube[9]);
  //}

  ///////////////end of hacky glitch thing for
  /// Vela///////////////////////////////////////////////////////

  // Convert priors to physical units (here only for timing parameters and
  // jumps)
  for (int p = 0; p < ((MNStruct *)globalcontext)->numFitTiming +
                          ((MNStruct *)globalcontext)->numFitJumps;
       p++) {
    if (((MNStruct *)globalcontext)->Dpriors[p][1] !=
        ((MNStruct *)globalcontext)->Dpriors[p][0]) {
      double val = 0;
      if ((((MNStruct *)globalcontext)->LDpriors[p][3]) == 0) {
        val = Cube[fitcount];
      }
      if ((((MNStruct *)globalcontext)->LDpriors[p][3]) == 1) {
        val = pow(10.0, Cube[fitcount]);
      }
      if ((((MNStruct *)globalcontext)->LDpriors[p][3]) == 2) {
        val = pow(10.0, Cube[fitcount]);
        uniformpriorterm += log(val);
      }

      LDparams[p] = val * (((MNStruct *)globalcontext)->LDpriors[p][1]) +
                    (((MNStruct *)globalcontext)->LDpriors[p][0]);

      if (((MNStruct *)globalcontext)->TempoFitNums[p][0] == param_sini &&
          ((MNStruct *)globalcontext)->usecosiprior == 1) {
        val = Cube[fitcount];
        LDparams[p] = std::sqrt(1.0 - val * val);
      }

      fitcount++;

    } else if (((MNStruct *)globalcontext)->Dpriors[p][1] ==
               ((MNStruct *)globalcontext)->Dpriors[p][0]) {
      LDparams[p] = ((MNStruct *)globalcontext)->Dpriors[p][0] *
                        (((MNStruct *)globalcontext)->LDpriors[p][1]) +
                    (((MNStruct *)globalcontext)->LDpriors[p][0]);
    }
  }

  pcount = 0;
  double phase = (double)LDparams[0];
  pcount++;
  for (int p = 1; p < ((MNStruct *)globalcontext)->numFitTiming; p++) {
    ((MNStruct *)globalcontext)
        ->pulse->param[((MNStruct *)globalcontext)->TempoFitNums[p][0]]
        .val[((MNStruct *)globalcontext)->TempoFitNums[p][1]] =
        LDparams[pcount];
    pcount++;
  }
  for (int p = 0; p < ((MNStruct *)globalcontext)->numFitJumps; p++) {
    ((MNStruct *)globalcontext)
        ->pulse->jumpVal[((MNStruct *)globalcontext)->TempoJumpNums[p]] =
        LDparams[pcount];
    pcount++;
  }

  if (TimetoMargin != numfit) {
    fastformBatsAll(((MNStruct *)globalcontext)->pulse,
                    ((MNStruct *)globalcontext)
                        ->numberpulsars); /* Form Barycentric arrival times */
    formResiduals(((MNStruct *)globalcontext)->pulse,
                  ((MNStruct *)globalcontext)->numberpulsars,
                  1); /* Form residuals */
  }

  for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
    Resvec[o] =
        (double)((MNStruct *)globalcontext)->pulse->obsn[o].residual + phase;
    // printf("Res: %i %g %g %g %g %g %g %g\n", o, Resvec[o], (double) phase,
    // ((MNStruct *)globalcontext)->Dpriors[0][0], ((MNStruct
    // *)globalcontext)->Dpriors[0][1], (double)((MNStruct
    // *)globalcontext)->LDpriors[0][1] , (double)((MNStruct
    // *)globalcontext)->LDpriors[0][0], (double)((MNStruct
    // *)globalcontext)->pulse->obsn[o].residual);
  }

  pcount = fitcount;
  if (((MNStruct *)globalcontext)->incStep > 0) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->incStep; i++) {
      int GrouptoFit = 0;
      double GroupStartTime = 0;

      GrouptoFit = floor(Cube[pcount]);
      pcount++;

      double GLength =
          ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][1] -
          ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][0];
      GroupStartTime =
          ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][0] +
          Cube[pcount] * GLength;
      pcount++;

      double StepAmp = Cube[pcount];
      pcount++;

      // printf("Step details: Group %i S %g F %g SS %g A %g \n", GrouptoFit,
      // ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][0],((MNStruct
      // *)globalcontext)->GroupStartTimes[GrouptoFit][1],GroupStartTime,StepAmp);

      for (int o1 = 0; o1 < ((MNStruct *)globalcontext)->pulse->nobs; o1++) {
        if (((MNStruct *)globalcontext)->pulse->obsn[o1].sat > GroupStartTime &&
            ((MNStruct *)globalcontext)->GroupNoiseFlags[o1] == GrouptoFit) {
          Resvec[o1] += StepAmp;
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Subtract
  /// GLitches/////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < ((MNStruct *)globalcontext)->incGlitch; i++) {
    double GlitchMJD = Cube[pcount];
    pcount++;

    double *GlitchAmps = new double[3];
    if (((MNStruct *)globalcontext)->incGlitchTerms == 1) {
      GlitchAmps[0] = Cube[pcount];
      pcount++;
    } else if (((MNStruct *)globalcontext)->incGlitchTerms == 2) {
      GlitchAmps[0] = Cube[pcount];
      pcount++;
      GlitchAmps[1] = Cube[pcount];
      pcount++;
    } else if (((MNStruct *)globalcontext)->incGlitchTerms == 3) {
      GlitchAmps[0] = Cube[pcount];
      pcount++;
      GlitchAmps[1] = pow(10.0, Cube[pcount]);  // Decay Amp
      pcount++;
      GlitchAmps[2] = pow(10.0, Cube[pcount]);  // Decay Timescale
      pcount++;
    }

    for (int o1 = 0; o1 < ((MNStruct *)globalcontext)->pulse->nobs; o1++) {
      if (((MNStruct *)globalcontext)->pulse->obsn[o1].bat > GlitchMJD) {
        if (((MNStruct *)globalcontext)->incGlitchTerms == 1) {
          long double arg = 0;
          arg =
              ((((MNStruct *)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
               ((MNStruct *)globalcontext)->pulse->param[param_f].val[0]) *
              86400.0;
          double darg = (double)arg;
          Resvec[o1] += GlitchAmps[0] * darg;

        } else if (((MNStruct *)globalcontext)->incGlitchTerms == 2) {
          for (int j = 0; j < ((MNStruct *)globalcontext)->incGlitchTerms;
               j++) {
            long double arg = 0;
            if (j == 0) {
              arg =
                  ((((MNStruct *)globalcontext)->pulse->obsn[o1].bat -
                    GlitchMJD) /
                   ((MNStruct *)globalcontext)->pulse->param[param_f].val[0]) *
                  86400.0;
            }
            if (j == 1) {
              arg = 0.5 *
                    pow((((MNStruct *)globalcontext)->pulse->obsn[o1].bat -
                         GlitchMJD) *
                            86400.0,
                        2) /
                    ((MNStruct *)globalcontext)->pulse->param[param_f].val[0];
            }
            double darg = (double)arg;
            Resvec[o1] += GlitchAmps[j] * darg;
          }
        } else if (((MNStruct *)globalcontext)->incGlitchTerms == 3) {
          long double arg = 0;
          double time =
              (((MNStruct *)globalcontext)->pulse->obsn[o1].bat - GlitchMJD);
          arg =
              ((((MNStruct *)globalcontext)->pulse->obsn[o1].bat - GlitchMJD) /
               ((MNStruct *)globalcontext)->pulse->param[param_f].val[0]) *
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
  /////////////////////////Get White Noise
  /// vector///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->numFitEFAC == 0) {
    EFAC = new double *[((MNStruct *)globalcontext)->EPolTerms];
    for (int n = 1; n <= ((MNStruct *)globalcontext)->EPolTerms; n++) {
      EFAC[n - 1] = new double[((MNStruct *)globalcontext)->systemcount];
      if (n == 1) {
        for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
          EFAC[n - 1][o] = 1;
        }
      } else {
        for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
          EFAC[n - 1][o] = 0;
        }
      }
    }
  } else if (((MNStruct *)globalcontext)->numFitEFAC == 1) {
    EFAC = new double *[((MNStruct *)globalcontext)->EPolTerms];
    for (int n = 1; n <= ((MNStruct *)globalcontext)->EPolTerms; n++) {
      EFAC[n - 1] = new double[((MNStruct *)globalcontext)->systemcount];
      if (n == 1) {
        for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
          EFAC[n - 1][o] = pow(10.0, Cube[pcount]);
          if (((MNStruct *)globalcontext)->EFACPriorType == 1) {
            uniformpriorterm += log(EFAC[n - 1][o]);
          }
        }
        pcount++;
      } else {
        for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
          EFAC[n - 1][o] = pow(10.0, Cube[pcount]);
        }
        pcount++;
      }
    }

  } else if (((MNStruct *)globalcontext)->numFitEFAC > 1) {
    EFAC = new double *[((MNStruct *)globalcontext)->EPolTerms];
    for (int n = 1; n <= ((MNStruct *)globalcontext)->EPolTerms; n++) {
      EFAC[n - 1] = new double[((MNStruct *)globalcontext)->systemcount];
      if (n == 1) {
        for (int p = 0; p < ((MNStruct *)globalcontext)->systemcount; p++) {
          EFAC[n - 1][p] = pow(10.0, Cube[pcount]);
          if (((MNStruct *)globalcontext)->EFACPriorType == 1) {
            uniformpriorterm += log(EFAC[n - 1][p]);
          }
          pcount++;
        }
      } else {
        for (int p = 0; p < ((MNStruct *)globalcontext)->systemcount; p++) {
          EFAC[n - 1][p] = pow(10.0, Cube[pcount]);
          pcount++;
        }
      }
    }
  }

  // printf("Equad %i \n", ((MNStruct *)globalcontext)->numFitEQUAD);
  if (((MNStruct *)globalcontext)->numFitEQUAD == 0) {
    EQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      EQUAD[o] = 0;
    }
  } else if (((MNStruct *)globalcontext)->numFitEQUAD == 1) {
    EQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      EQUAD[o] = pow(10.0, 2 * Cube[pcount]);
      if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
        uniformpriorterm += log(pow(10.0, Cube[pcount]));
      }
    }
    pcount++;
  } else if (((MNStruct *)globalcontext)->numFitEQUAD > 1) {
    EQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      if (((MNStruct *)globalcontext)->includeEQsys[o] == 1) {
        // printf("Cube: %i %i %g \n", o, pcount, Cube[pcount]);
        EQUAD[o] = pow(10.0, 2 * Cube[pcount]);
        if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
          uniformpriorterm += log(pow(10.0, Cube[pcount]));
        }
        pcount++;
      } else {
        EQUAD[o] = 0;
      }
      // printf("Equad? %i %g \n", o, EQUAD[o]);
    }
  }

  double *SQUAD;
  if (((MNStruct *)globalcontext)->incShannonJitter == 0) {
    SQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      SQUAD[o] = 0;
    }
  } else if (((MNStruct *)globalcontext)->incShannonJitter == 1) {
    SQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      SQUAD[o] = pow(10.0, 2 * Cube[pcount]);
      if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
        uniformpriorterm += log(pow(10.0, Cube[pcount]));
      }
    }
    pcount++;
  } else if (((MNStruct *)globalcontext)->incShannonJitter > 1) {
    SQUAD = new double[((MNStruct *)globalcontext)->systemcount];
    for (int o = 0; o < ((MNStruct *)globalcontext)->systemcount; o++) {
      SQUAD[o] = pow(10.0, 2 * Cube[pcount]);
      if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
        uniformpriorterm += log(pow(10.0, Cube[pcount]));
      }

      pcount++;
    }
  }

  double *ECORRPrior;
  if (((MNStruct *)globalcontext)->incNGJitter > 0) {
    double *ECorrCoeffs = new double[((MNStruct *)globalcontext)->incNGJitter];
    for (int i = 0; i < ((MNStruct *)globalcontext)->incNGJitter; i++) {
      ECorrCoeffs[i] = pow(10.0, 2 * Cube[pcount]);
      if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
        uniformpriorterm += log(pow(10.0, Cube[pcount]));
      }
      pcount++;
    }
    ECORRPrior = new double[((MNStruct *)globalcontext)->numNGJitterEpochs];
    for (int i = 0; i < ((MNStruct *)globalcontext)->numNGJitterEpochs; i++) {
      ECORRPrior[i] =
          ECorrCoeffs[((MNStruct *)globalcontext)->NGJitterSysFlags[i]];
    }

    delete[] ECorrCoeffs;
  }

  double *SECORRPrior;
  if (((MNStruct *)globalcontext)->incNGSJitter > 0) {
    double *SECorrCoeffs =
        new double[((MNStruct *)globalcontext)->incNGSJitter];
    for (int i = 0; i < ((MNStruct *)globalcontext)->incNGSJitter; i++) {
      SECorrCoeffs[i] = pow(10.0, 2 * Cube[pcount]);
      if (((MNStruct *)globalcontext)->EQUADPriorType == 1) {
        uniformpriorterm += log(pow(10.0, Cube[pcount]));
      }
      pcount++;
    }
    SECORRPrior = new double[((MNStruct *)globalcontext)->numNGSJitterEpochs];
    for (int i = 0; i < ((MNStruct *)globalcontext)->numNGSJitterEpochs; i++) {
      SECORRPrior[i] =
          SECorrCoeffs[((MNStruct *)globalcontext)->NGSJitterSysFlags[i]];
    }

    delete[] SECorrCoeffs;
  }

  double DMEQUAD = 0;
  if (((MNStruct *)globalcontext)->incDMEQUAD > 0) {
    DMEQUAD = pow(10.0, Cube[pcount]);
    pcount++;
  }

  double SolarWind = 0;
  double WhiteSolarWind = 0;

  if (((MNStruct *)globalcontext)->FitSolarWind == 1) {
    SolarWind = Cube[pcount];
    pcount++;

    // printf("Solar Wind: %g \n", SolarWind);
  }
  if (((MNStruct *)globalcontext)->FitWhiteSolarWind == 1) {
    WhiteSolarWind = pow(10.0, Cube[pcount]);
    pcount++;

    // printf("White Solar Wind: %g \n", WhiteSolarWind);
  }

  double *Noise;
  double *BATvec;
  Noise = new double[((MNStruct *)globalcontext)->pulse->nobs];
  // BATvec=new double[((MNStruct *)globalcontext)->pulse->nobs];

  // for(int o=0;o<((MNStruct *)globalcontext)->pulse->nobs; o++){
  // BATvec[o]=(double)((MNStruct *)globalcontext)->pulse->obsn[o].bat;
  //}

  double DMKappa = 2.410 * pow(10.0, -16);
  if (((MNStruct *)globalcontext)->whitemodel == 0) {
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      double EFACterm = 0;
      double noiseval = 0;
      double ShannonJitterTerm = 0;

      double SWTerm = WhiteSolarWind *
                      ((MNStruct *)globalcontext)->pulse->obsn[o].tdis2 /
                      ((MNStruct *)globalcontext)->pulse->ne_sw;
      double DMEQUADTerm =
          DMEQUAD /
          (DMKappa *
           pow((double)((MNStruct *)globalcontext)->pulse->obsn[o].freqSSB, 2));

      if (((MNStruct *)globalcontext)->useOriginalErrors == 0) {
        noiseval = ((MNStruct *)globalcontext)->pulse->obsn[o].toaErr;
      } else if (((MNStruct *)globalcontext)->useOriginalErrors == 1) {
        noiseval = ((MNStruct *)globalcontext)->pulse->obsn[o].origErr;
      }

      for (int n = 1; n <= ((MNStruct *)globalcontext)->EPolTerms; n++) {
        EFACterm =
            EFACterm +
            pow((noiseval * pow(10.0, -6)) / pow(pow(10.0, -7), n - 1), n) *
                EFAC[n - 1][((MNStruct *)globalcontext)->sysFlags[o]];
      }

      if (((MNStruct *)globalcontext)->incShannonJitter > 0) {
        ShannonJitterTerm = SQUAD[((MNStruct *)globalcontext)->sysFlags[o]] *
                            ((MNStruct *)globalcontext)->pulse->obsn[o].tobs /
                            1000.0;
      }
      // printf("Noise: %i %g %g %g %g %g \n", EFACterm, EQUAD[((MNStruct
      // *)globalcontext)->sysFlags[o]], ShannonJitterTerm, SWTerm,
      // DMEQUADTerm);
      Noise[o] =
          1.0 /
          (pow(EFACterm, 2) + EQUAD[((MNStruct *)globalcontext)->sysFlags[o]] +
           ShannonJitterTerm + SWTerm * SWTerm + DMEQUADTerm * DMEQUADTerm);
    }

  } else if (((MNStruct *)globalcontext)->whitemodel == 1) {
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      Noise[o] =
          1.0 / (EFAC[0][((MNStruct *)globalcontext)->sysFlags[o]] *
                 EFAC[0][((MNStruct *)globalcontext)->sysFlags[o]] *
                 (pow(((((MNStruct *)globalcontext)->pulse->obsn[o].toaErr) *
                       pow(10.0, -6)),
                      2) +
                  EQUAD[((MNStruct *)globalcontext)->sysFlags[o]]));
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Initialise
  /// TotalMatrix////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  int totalsize = ((MNStruct *)globalcontext)->totalsize;

  double *TotalMatrix = ((MNStruct *)globalcontext)->StoredTMatrix;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Form the Design
  /// Matrix////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  if (TimetoMargin != ((MNStruct *)globalcontext)->numFitTiming +
                          ((MNStruct *)globalcontext)->numFitJumps) {
    getCustomDVectorLike(globalcontext, TotalMatrix,
                         ((MNStruct *)globalcontext)->pulse->nobs, TimetoMargin,
                         totalsize);
    vector_dgesvd(TotalMatrix, ((MNStruct *)globalcontext)->pulse->nobs,
                  TimetoMargin);
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////Set up
  /// Coefficients///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  double maxtspan = ((MNStruct *)globalcontext)->Tspan;
  double averageTSamp = 2 * maxtspan / ((MNStruct *)globalcontext)->pulse->nobs;

  int FitRedCoeff = 2 * (((MNStruct *)globalcontext)->numFitRedCoeff);
  int FitDMCoeff = 2 * (((MNStruct *)globalcontext)->numFitDMCoeff);
  int FitScatCoeff = 2 * (((MNStruct *)globalcontext)->numFitScatCoeff);
  int FitBandCoeff = 2 * (((MNStruct *)globalcontext)->numFitBandNoiseCoeff);
  int FitGroupNoiseCoeff =
      2 * ((MNStruct *)globalcontext)->numFitGroupNoiseCoeff;

  int totCoeff = ((MNStruct *)globalcontext)->totCoeff;

  double *powercoeff = new double[totCoeff];
  for (int o = 0; o < totCoeff; o++) {
    powercoeff[o] = 0;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Red
  /// Noise///////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double *freqs = new double[totCoeff];
  double *DMVec = new double[((MNStruct *)globalcontext)->pulse->nobs];

  double freqdet = 0;
  int startpos = 0;

  if (((MNStruct *)globalcontext)->incRED > 0 ||
      ((MNStruct *)globalcontext)->incGWB == 1) {
    if (((MNStruct *)globalcontext)->FitLowFreqCutoff == 1) {
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
        ((MNStruct *)globalcontext)->sampleFreq[freqStartpoint] =
            pow(10.0, sampledFLow + i * RedLogSampledDiff / LogLowFreqs);
        freqStartpoint++;
      }

      for (int i = 0; i < FitRedCoeff / 2 - LogLowFreqs; i++) {
        ((MNStruct *)globalcontext)->sampleFreq[freqStartpoint] =
            i + RedMidFreq;
        freqStartpoint++;
      }
    }

    if (((MNStruct *)globalcontext)->FitLowFreqCutoff == 2) {
      double fLow = pow(10.0, Cube[pcount]);
      pcount++;

      for (int i = 0; i < FitRedCoeff / 2; i++) {
        ((MNStruct *)globalcontext)->sampleFreq[i] = ((double)(i + 1)) * fLow;
      }
    }

    for (int i = 0;
         i < FitRedCoeff / 2 - ((MNStruct *)globalcontext)->incFloatRed; i++) {
      freqs[startpos + i] =
          (double)((MNStruct *)globalcontext)->sampleFreq[i] / maxtspan;
      freqs[startpos + i + FitRedCoeff / 2] = freqs[startpos + i];
    }
  }

  if (((MNStruct *)globalcontext)->storeFMatrices == 0) {
    for (int i = 0; i < FitRedCoeff / 2; i++) {
      for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
        double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

        TotalMatrix[k + (i + TimetoMargin + startpos) *
                            ((MNStruct *)globalcontext)->pulse->nobs] =
            cos(2 * M_PI * freqs[i] * time);
        TotalMatrix[k + (i + FitRedCoeff / 2 + TimetoMargin + startpos) *
                            ((MNStruct *)globalcontext)->pulse->nobs] =
            sin(2 * M_PI * freqs[i] * time);
      }
    }
  }

  if (((MNStruct *)globalcontext)->incRED == 2) {
    for (int i = 0; i < FitRedCoeff / 2; i++) {
      int pnum = pcount;
      double pc = Cube[pcount];

      if (((MNStruct *)globalcontext)->RedPriorType == 1) {
        uniformpriorterm += log(pow(10.0, pc));
      }

      powercoeff[i] = pow(10.0, 2 * pc);
      powercoeff[i + FitRedCoeff / 2] = powercoeff[i];
      pcount++;
    }

    startpos = FitRedCoeff;

  } else if (((MNStruct *)globalcontext)->incRED == 3 ||
             ((MNStruct *)globalcontext)->incRED == 4) {
    for (int pl = 0; pl < ((MNStruct *)globalcontext)->numFitRedPL; pl++) {
      double Tspan = maxtspan;
      double f1yr = 1.0 / 3.16e7;

      if (((MNStruct *)globalcontext)->FitLowFreqCutoff == 2) {
        Tspan = Tspan / ((MNStruct *)globalcontext)->sampleFreq[0];
      }

      double redamp = Cube[pcount];
      pcount++;
      double redindex = Cube[pcount];
      pcount++;

      double cornerfreq = 0;
      if (((MNStruct *)globalcontext)->incRED == 4) {
        cornerfreq = pow(10.0, Cube[pcount]) / Tspan;
        pcount++;
      }

      redamp = pow(10.0, redamp);
      if (((MNStruct *)globalcontext)->RedPriorType == 1) {
        uniformpriorterm += log(redamp);
      }

      double Agw = redamp;
      for (int i = 0;
           i < FitRedCoeff / 2 - ((MNStruct *)globalcontext)->incFloatRed;
           i++) {
        double rho = 0;
        if (((MNStruct *)globalcontext)->incRED == 3) {
          rho = (Agw * Agw / 12.0 / (M_PI * M_PI)) * pow(f1yr, (-3)) *
                pow(freqs[i] * 365.25, (-redindex)) / (Tspan * 24 * 60 * 60);
        }
        if (((MNStruct *)globalcontext)->incRED == 4) {
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
    if (((MNStruct *)globalcontext)->varyRedCoeff == 1) {
      coefftovary = int(pow(10.0, Cube[pcount])) - 1;
      pcount++;
      amptovary = pow(10.0, Cube[pcount]);
      pcount++;

      powercoeff[coefftovary] = amptovary;
      powercoeff[coefftovary + FitRedCoeff / 2] = amptovary;
    }

    startpos = FitRedCoeff;
  }

  if (((MNStruct *)globalcontext)->incGWB == 1) {
    double GWBAmp = pow(10.0, Cube[pcount]);
    pcount++;
    uniformpriorterm += log(GWBAmp);
    double Tspan = maxtspan;

    if (((MNStruct *)globalcontext)->FitLowFreqCutoff == 2) {
      Tspan = Tspan / ((MNStruct *)globalcontext)->sampleFreq[0];
    }

    double f1yr = 1.0 / 3.16e7;
    for (int i = 0;
         i < FitRedCoeff / 2 - ((MNStruct *)globalcontext)->incFloatRed; i++) {
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
  /////////////////////////Subtract  Red Shape
  /// Events////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  int MarginRedShapelets = ((MNStruct *)globalcontext)->MarginRedShapeCoeff;
  int totalredshapecoeff = 0;
  int badshape = 0;

  double **RedShapeletMatrix;

  double globalwidth = 0;
  if (((MNStruct *)globalcontext)->incRedShapeEvent != 0) {
    if (MarginRedShapelets == 1) {
      badshape = 1;

      totalredshapecoeff = ((MNStruct *)globalcontext)->numRedShapeCoeff *
                           ((MNStruct *)globalcontext)->incRedShapeEvent;

      RedShapeletMatrix =
          new double *[((MNStruct *)globalcontext)->pulse->nobs];
      for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
        RedShapeletMatrix[i] = new double[totalredshapecoeff];
        for (int j = 0; j < totalredshapecoeff; j++) {
          RedShapeletMatrix[i][j] = 0;
        }
      }
    }

    for (int i = 0; i < ((MNStruct *)globalcontext)->incRedShapeEvent; i++) {
      int numRedShapeCoeff = ((MNStruct *)globalcontext)->numRedShapeCoeff;

      double EventPos = Cube[pcount];
      pcount++;
      double EventWidth = Cube[pcount];
      pcount++;

      globalwidth = EventWidth;

      double *RedshapeVec = new double[numRedShapeCoeff];

      double *RedshapeNorm = new double[numRedShapeCoeff];
      for (int c = 0; c < numRedShapeCoeff; c++) {
        RedshapeNorm[c] = 1.0 / std::sqrt(std::sqrt(2.0 * M_PI) * pow(2.0, c) *
                                          iter_factorial(c));
      }

      for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
        double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

        double HVal = (time - EventPos) / (std::sqrt(2.0) * EventWidth);
        TNothpl(numRedShapeCoeff, HVal, RedshapeVec);

        for (int c = 0; c < numRedShapeCoeff; c++) {
          RedShapeletMatrix[k][i * numRedShapeCoeff + c] =
              RedshapeNorm[c] * RedshapeVec[c] *
              exp(-0.5 * pow((time - EventPos) / EventWidth, 2));
          if (RedShapeletMatrix[k][i * numRedShapeCoeff + c] > 0.01) {
            badshape = 0;
          }
        }
      }

      delete[] RedshapeVec;
      delete[] RedshapeNorm;
    }

    if (MarginRedShapelets == 0) {
      for (int i = 0; i < ((MNStruct *)globalcontext)->incRedShapeEvent; i++) {
        int numRedShapeCoeff = ((MNStruct *)globalcontext)->numRedShapeCoeff;

        double EventPos = Cube[pcount];
        pcount++;
        double EventWidth = Cube[pcount];
        pcount++;

        globalwidth = EventWidth;

        double *Redshapecoeff = new double[numRedShapeCoeff];
        double *RedshapeVec = new double[numRedShapeCoeff];
        for (int c = 0; c < numRedShapeCoeff; c++) {
          Redshapecoeff[c] = Cube[pcount];
          pcount++;
        }

        double *RedshapeNorm = new double[numRedShapeCoeff];
        for (int c = 0; c < numRedShapeCoeff; c++) {
          RedshapeNorm[c] = 1.0 / std::sqrt(std::sqrt(2.0 * M_PI) *
                                            pow(2.0, c) * iter_factorial(c));
        }

        for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
          double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

          double HVal = (time - EventPos) / (std::sqrt(2.0) * EventWidth);
          TNothpl(numRedShapeCoeff, HVal, RedshapeVec);
          double Redsignal = 0;
          for (int c = 0; c < numRedShapeCoeff; c++) {
            Redsignal += RedshapeNorm[c] * RedshapeVec[c] * Redshapecoeff[c];
          }

          Resvec[k] -=
              Redsignal * exp(-0.5 * pow((time - EventPos) / EventWidth, 2));

          // printf("Shape Sig: %i %.10Lg %g \n", k, ((MNStruct
          // *)globalcontext)->pulse->obsn[k].bat,
          // Redsignal*exp(-0.5*pow((time-EventPos)/EventWidth, 2)));
        }

        delete[] Redshapecoeff;
        delete[] RedshapeVec;
        delete[] RedshapeNorm;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Subtract Sine
  /// Wave///////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incsinusoid == 1) {
    double sineamp = pow(10.0, Cube[pcount]);
    pcount++;
    double sinephase = Cube[pcount];
    pcount++;
    double sinefreq = pow(10.0, Cube[pcount]) / maxtspan;
    pcount++;
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      Resvec[o] -=
          sineamp *
          sin(2 * M_PI * sinefreq *
                  (double)((MNStruct *)globalcontext)->pulse->obsn[o].bat +
              sinephase);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////DM
  /// Variations////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incDM > 0) {
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      DMVec[o] =
          1.0 /
          (DMKappa *
           pow((double)((MNStruct *)globalcontext)->pulse->obsn[o].freqSSB, 2));
    }

    for (int i = 0; i < FitDMCoeff / 2; i++) {
      freqs[startpos + i] =
          ((MNStruct *)globalcontext)
              ->sampleFreq[startpos / 2 -
                           ((MNStruct *)globalcontext)->incFloatRed + i] /
          maxtspan;
      freqs[startpos + i + FitDMCoeff / 2] = freqs[startpos + i];

      if (((MNStruct *)globalcontext)->storeFMatrices == 0) {
        for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
          double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

          TotalMatrix[k + (i + TimetoMargin + startpos) *
                              ((MNStruct *)globalcontext)->pulse->nobs] =
              cos(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
          TotalMatrix[k + (i + FitDMCoeff / 2 + TimetoMargin + startpos) *
                              ((MNStruct *)globalcontext)->pulse->nobs] =
              sin(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
        }
      }
    }
  }

  if (((MNStruct *)globalcontext)->incDM == 2) {
    for (int i = 0; i < FitDMCoeff / 2; i++) {
      int pnum = pcount;
      double pc = Cube[pcount];

      powercoeff[startpos + i] = pow(10.0, 2 * pc);
      powercoeff[startpos + i + FitDMCoeff / 2] = powercoeff[startpos + i];
      freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
      pcount++;
    }
    startpos += FitDMCoeff;

  } else if (((MNStruct *)globalcontext)->incDM == 3) {
    for (int pl = 0; pl < ((MNStruct *)globalcontext)->numFitDMPL; pl++) {
      double DMamp = Cube[pcount];
      pcount++;
      double DMindex = Cube[pcount];
      pcount++;

      double Tspan = maxtspan;
      double f1yr = 1.0 / 3.16e7;

      DMamp = pow(10.0, DMamp);
      if (((MNStruct *)globalcontext)->DMPriorType == 1) {
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
    if (((MNStruct *)globalcontext)->varyDMCoeff == 1) {
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
  /////////////////////////Scattering variations
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  double *ScatVec = new double[((MNStruct *)globalcontext)->pulse->nobs];

  if (((MNStruct *)globalcontext)->incScat > 0) {
    if (((MNStruct *)globalcontext)->storeFMatrices == 0)
      for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++)
        ScatVec[o] =
            1. /
            (pow(
                (double)((MNStruct *)globalcontext)->pulse->obsn[o].freqSSB /
                    1400e6,
                4));  // Referenced to 1.4 Ghz to be consistent with Enterprise;
                      // Invert to be consistent with Tempo2 - 20220810 GD/AP

    for (int i = 0; i < FitScatCoeff / 2; i++) {
      freqs[startpos + i] =
          ((MNStruct *)globalcontext)
              ->sampleFreq[startpos / 2 -
                           ((MNStruct *)globalcontext)->incFloatRed + i] /
          maxtspan;
      freqs[startpos + i + FitScatCoeff / 2] = freqs[startpos + i];

      if (((MNStruct *)globalcontext)->storeFMatrices == 0) {
        for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
          double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

          TotalMatrix[k + (i + TimetoMargin + startpos) *
                              ((MNStruct *)globalcontext)->pulse->nobs] =
              cos(2 * M_PI * freqs[startpos + i] * time) * ScatVec[k];
          TotalMatrix[k + (i + FitScatCoeff / 2 + TimetoMargin + startpos) *
                              ((MNStruct *)globalcontext)->pulse->nobs] =
              sin(2 * M_PI * freqs[startpos + i] * time) * ScatVec[k];
        }
      }
    }

    for (int pl = 0; pl < ((MNStruct *)globalcontext)->numFitScatPL; pl++) {
      double ScatAmp = Cube[pcount];
      pcount++;
      double ScatSlope = Cube[pcount];
      pcount++;

      double Tspan = maxtspan;
      double f1yr = 1.0 / 3.16e7;

      ScatAmp = pow(10.0, ScatAmp);
      if (((MNStruct *)globalcontext)->ScatPriorType == 1) {
        uniformpriorterm += log(ScatAmp);
      }
      for (int i = 0; i < FitScatCoeff / 2; i++) {
        double rho = (ScatAmp * ScatAmp) / 12. / M_PI / M_PI * pow(f1yr, (-3)) *
                     pow(freqs[startpos + i] * 365.25, (-ScatSlope)) /
                     (maxtspan * 24 * 60 * 60);
        powercoeff[startpos + i] += rho;
        powercoeff[startpos + i + FitScatCoeff / 2] += rho;
      }
    }

    for (int i = 0; i < FitScatCoeff / 2; i++) {
      freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
    }
    startpos += FitScatCoeff;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Subtract  DM Shape
  /// Events////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incDMShapeEvent != 0) {
    if (((MNStruct *)globalcontext)->incDM == 0) {
      for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
        DMVec[o] =
            1.0 /
            (DMKappa *
             pow((double)((MNStruct *)globalcontext)->pulse->obsn[o].freqSSB,
                 2));
      }
    }
    for (int i = 0; i < ((MNStruct *)globalcontext)->incDMShapeEvent; i++) {
      int numDMShapeCoeff = ((MNStruct *)globalcontext)->numDMShapeCoeff;

      double EventPos = Cube[pcount];
      pcount++;
      double EventWidth = Cube[pcount];
      pcount++;

      globalwidth = EventWidth;

      double *DMshapecoeff = new double[numDMShapeCoeff];
      double *DMshapeVec = new double[numDMShapeCoeff];
      for (int c = 0; c < numDMShapeCoeff; c++) {
        DMshapecoeff[c] = Cube[pcount];
        pcount++;
      }

      double *DMshapeNorm = new double[numDMShapeCoeff];
      for (int c = 0; c < numDMShapeCoeff; c++) {
        DMshapeNorm[c] = 1.0 / std::sqrt(std::sqrt(2.0 * M_PI) * pow(2.0, c) *
                                         iter_factorial(c));
      }

      for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
        double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

        double HVal = (time - EventPos) / (std::sqrt(2.0) * EventWidth);
        TNothpl(numDMShapeCoeff, HVal, DMshapeVec);
        double DMsignal = 0;
        for (int c = 0; c < numDMShapeCoeff; c++) {
          DMsignal +=
              DMshapeNorm[c] * DMshapeVec[c] * DMshapecoeff[c] * DMVec[k];
        }

        Resvec[k] -=
            DMsignal * exp(-0.5 * pow((time - EventPos) / EventWidth, 2));
      }

      delete[] DMshapecoeff;
      delete[] DMshapeVec;
      delete[] DMshapeNorm;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Subtract  DM Scatter Shape
  /// Events///////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incDMScatterShapeEvent != 0) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->incDMScatterShapeEvent;
         i++) {
      int numDMShapeCoeff = ((MNStruct *)globalcontext)->numDMScatterShapeCoeff;

      double EventPos = Cube[pcount];
      pcount++;
      double EventWidth = globalwidth;  // Cube[pcount];
      pcount++;
      double EventFScale = Cube[pcount];
      pcount++;

      // EventWidth=globalwidth;

      double *DMshapecoeff = new double[numDMShapeCoeff];
      double *DMshapeVec = new double[numDMShapeCoeff];
      for (int c = 0; c < numDMShapeCoeff; c++) {
        DMshapecoeff[c] = Cube[pcount];
        pcount++;
      }

      double *DMshapeNorm = new double[numDMShapeCoeff];
      for (int c = 0; c < numDMShapeCoeff; c++) {
        DMshapeNorm[c] = 1.0 / std::sqrt(std::sqrt(2.0 * M_PI) * pow(2.0, c) *
                                         iter_factorial(c));
      }

      for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
        double time = (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;

        double HVal = (time - EventPos) / (std::sqrt(2.0) * EventWidth);
        TNothpl(numDMShapeCoeff, HVal, DMshapeVec);
        double DMsignal = 0;
        for (int c = 0; c < numDMShapeCoeff; c++) {
          DMsignal += DMshapeNorm[c] * DMshapeVec[c] * DMshapecoeff[c] *
                      pow(DMVec[k], EventFScale / 2.0);
        }

        Resvec[k] -=
            DMsignal * exp(-0.5 * pow((time - EventPos) / EventWidth, 2));
      }

      delete[] DMshapecoeff;
      delete[] DMshapeVec;
      delete[] DMshapeNorm;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Subtract Yearly
  /// DM//////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->yearlyDM == 1) {
    double yearlyamp = pow(10.0, Cube[pcount]);
    pcount++;
    double yearlyphase = Cube[pcount];
    pcount++;
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      Resvec[o] -=
          yearlyamp *
          sin((2 * M_PI / 365.25) *
                  (double)((MNStruct *)globalcontext)->pulse->obsn[o].bat +
              yearlyphase) *
          DMVec[o];
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Subtract Solar
  /// Wind//////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->FitSolarWind == 1) {
    for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
      Resvec[o] -= (SolarWind - ((MNStruct *)globalcontext)->pulse->ne_sw) *
                   ((MNStruct *)globalcontext)->pulse->obsn[o].tdis2;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Band
  /// DM/////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incBandNoise > 0) {
    if (((MNStruct *)globalcontext)->incDM == 0) {
      for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
        DMVec[o] =
            1.0 /
            (DMKappa *
             pow((double)((MNStruct *)globalcontext)->pulse->obsn[o].freqSSB,
                 2));
      }
    }

    for (int b = 0; b < ((MNStruct *)globalcontext)->incBandNoise; b++) {
      double startfreq = ((MNStruct *)globalcontext)->FitForBand[b][0];
      double stopfreq = ((MNStruct *)globalcontext)->FitForBand[b][1];
      double BandScale = ((MNStruct *)globalcontext)->FitForBand[b][2];
      int BandPriorType = ((MNStruct *)globalcontext)->FitForBand[b][3];

      double Bandamp = Cube[pcount];
      pcount++;
      double Bandindex = Cube[pcount];
      pcount++;

      double Tspan = maxtspan;
      double f1yr = 1.0 / 3.16e7;

      Bandamp = pow(10.0, Bandamp);
      if (BandPriorType == 1) {
        uniformpriorterm += log(Bandamp);
      }
      for (int i = 0; i < FitBandCoeff / 2; i++) {
        freqs[startpos + i] = ((double)(i + 1)) / maxtspan;
        freqs[startpos + i + FitBandCoeff / 2] = freqs[startpos + i];

        double rho = (Bandamp * Bandamp / 12.0 / (M_PI * M_PI)) *
                     pow(f1yr, (-3)) *
                     pow(freqs[startpos + i] * 365.25, (-Bandindex)) /
                     (maxtspan * 24 * 60 * 60);
        powercoeff[startpos + i] += rho;
        powercoeff[startpos + i + FitBandCoeff / 2] += rho;
      }

      for (int i = 0; i < FitBandCoeff / 2; i++) {
        freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
      }
      if (((MNStruct *)globalcontext)->storeFMatrices == 0) {
        for (int i = 0; i < FitBandCoeff / 2; i++) {
          for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
            if (((MNStruct *)globalcontext)->pulse->obsn[k].freq > startfreq &&
                ((MNStruct *)globalcontext)->pulse->obsn[k].freq < stopfreq) {
              double time =
                  (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;
              TotalMatrix[k + (i + TimetoMargin + startpos) *
                                  ((MNStruct *)globalcontext)->pulse->nobs] =
                  cos(2 * M_PI * freqs[startpos + i] * time);
              TotalMatrix[k + (i + TimetoMargin + startpos + FitBandCoeff / 2) *
                                  ((MNStruct *)globalcontext)->pulse->nobs] =
                  sin(2 * M_PI * freqs[startpos + i] * time);
            } else {
              TotalMatrix[k + (i + TimetoMargin + startpos) *
                                  ((MNStruct *)globalcontext)->pulse->nobs] = 0;
              TotalMatrix[k + (i + TimetoMargin + startpos + FitBandCoeff / 2) *
                                  ((MNStruct *)globalcontext)->pulse->nobs] = 0;
            }
          }
        }
      }

      startpos += FitBandCoeff;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Add Group
  /// Noise/////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incGroupNoise > 0) {
    for (int g = 0; g < ((MNStruct *)globalcontext)->incGroupNoise; g++) {
      int GrouptoFit = 0;
      double GroupStartTime = 0;
      double GroupStopTime = 0;
      double GroupTSpan = 0;
      if (((MNStruct *)globalcontext)->FitForGroup[g][0] == -1) {
        GrouptoFit = floor(Cube[pcount]);
        pcount++;
        // printf("Fit for group %i \n", GrouptoFit);
      } else {
        GrouptoFit = ((MNStruct *)globalcontext)->FitForGroup[g][0];
      }

      if (((MNStruct *)globalcontext)->FitForGroup[g][1] == 1) {
        double GLength =
            ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][1] -
            ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][0];
        GroupStartTime =
            ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][0] +
            Cube[pcount] * GLength;

        double LengthLeft =
            ((MNStruct *)globalcontext)->GroupStartTimes[GrouptoFit][1] -
            GroupStartTime;

        pcount++;
        GroupStopTime = GroupStartTime + Cube[pcount] * LengthLeft;
        pcount++;
        GroupTSpan = GroupStopTime - GroupStartTime;

        // printf("Start Stop %i %g %g %g %g \n", GrouptoFit, ((MNStruct
        // *)globalcontext)->GroupStartTimes[GrouptoFit][0] , ((MNStruct
        // *)globalcontext)->GroupStartTimes[GrouptoFit][1] , GLength,
        // GroupStartTime);

      } else {
        GroupStartTime = 0;
        GroupStopTime = 10000000.0;
        GroupTSpan = maxtspan;
      }

      double GroupScale = ((MNStruct *)globalcontext)->FitForGroup[g][5];

      double GroupAmp = Cube[pcount];
      pcount++;
      double GroupIndex = Cube[pcount];
      pcount++;

      double Tspan = maxtspan;
      double f1yr = 1.0 / 3.16e7;

      GroupAmp = pow(10.0, GroupAmp);
      if (((MNStruct *)globalcontext)->DMPriorType == 1) {
        uniformpriorterm += log(GroupAmp);
      }

      for (int i = 0; i < FitGroupNoiseCoeff / 2; i++) {
        freqs[startpos + i] = ((double)(i + 1.0)) / GroupTSpan;  // maxtspan;
        freqs[startpos + i + FitGroupNoiseCoeff / 2] = freqs[startpos + i];

        double rho = (GroupAmp * GroupAmp / 12.0 / (M_PI * M_PI)) *
                     pow(f1yr, (-3)) *
                     pow(freqs[startpos + i] * 365.25, (-GroupIndex)) /
                     (maxtspan * 24 * 60 * 60);
        powercoeff[startpos + i] += rho;
        powercoeff[startpos + i + FitGroupNoiseCoeff / 2] += rho;
      }

      for (int i = 0; i < FitGroupNoiseCoeff / 2; i++) {
        freqdet = freqdet + 2 * log(powercoeff[startpos + i]);
      }

      double sqrtDMKappa = std::sqrt(DMKappa);
      for (int i = 0; i < FitGroupNoiseCoeff / 2; i++) {
        for (int k = 0; k < ((MNStruct *)globalcontext)->pulse->nobs; k++) {
          if (((MNStruct *)globalcontext)->GroupNoiseFlags[k] == GrouptoFit &&
              (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat >
                  GroupStartTime &&
              (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat <
                  GroupStopTime) {
            //		if(((MNStruct *)globalcontext)->GroupNoiseFlags[k] ==
            // GrouptoFit){
            double time =
                (double)((MNStruct *)globalcontext)->pulse->obsn[k].bat;
            double Fscale =
                1.0 /
                (std::pow(sqrtDMKappa * ((double)((MNStruct *)globalcontext)
                                             ->pulse->obsn[k]
                                             .freqSSB),
                          GroupScale));
            // printf("Scale: %i %i %g %g %g \n", g, k, GroupScale, Fscale);
            TotalMatrix[k + (i + TimetoMargin + startpos) *
                                ((MNStruct *)globalcontext)->pulse->nobs] =
                cos(2 * M_PI * freqs[startpos + i] * time) * Fscale;
            TotalMatrix[k +
                        (i + TimetoMargin + startpos + FitGroupNoiseCoeff / 2) *
                            ((MNStruct *)globalcontext)->pulse->nobs] =
                sin(2 * M_PI * freqs[startpos + i] * time) * Fscale;
          } else {
            TotalMatrix[k + (i + TimetoMargin + startpos) *
                                ((MNStruct *)globalcontext)->pulse->nobs] = 0;
            TotalMatrix[k +
                        (i + TimetoMargin + startpos + FitGroupNoiseCoeff / 2) *
                            ((MNStruct *)globalcontext)->pulse->nobs] = 0;
          }
        }
      }

      startpos = startpos + FitGroupNoiseCoeff;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Add ECORR
  /// Coeffs////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (((MNStruct *)globalcontext)->incNGJitter > 0) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->numNGJitterEpochs; i++) {
      powercoeff[startpos + i] = ECORRPrior[i];
      freqdet = freqdet + log(ECORRPrior[i]);
    }
  }

  if (((MNStruct *)globalcontext)->incNGSJitter > 0) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->numNGSJitterEpochs; i++) {
      powercoeff[startpos + i] = SECORRPrior[i];
      freqdet = freqdet + log(SECORRPrior[i]);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Get Time domain
  /// likelihood//////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double tdet = 0;
  double timelike = 0;

  for (int o = 0; o < ((MNStruct *)globalcontext)->pulse->nobs; o++) {
    // printf("Res: %i  %g %g \n", o, Resvec[o], sqrt(Noise[o]));
    timelike += Resvec[o] * Resvec[o] * Noise[o];
    tdet -= log(Noise[o]);
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////Form Total
  /// Matrices////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
    for (int j = 0; j < totalredshapecoeff; j++) {
      TotalMatrix[i + (j + TimetoMargin + totCoeff) *
                          ((MNStruct *)globalcontext)->pulse->nobs] =
          RedShapeletMatrix[i][j];
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////Do
  /// Algebra/////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  logtchk("Starting algebra");
  int savememory = 0;

  double *NTd = new double[totalsize];
  double *TNT = new double[totalsize * totalsize];
  double *NT;
  if (savememory == 0) {
    NT = new double[((MNStruct *)globalcontext)->pulse->nobs * totalsize];
    std::memcpy(
        NT, TotalMatrix,
        ((MNStruct *)globalcontext)->pulse->nobs * totalsize * sizeof(double));

    for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
      for (int j = 0; j < totalsize; j++) {
        NT[i + j * ((MNStruct *)globalcontext)->pulse->nobs] *= Noise[i];
      }
    }

    vector_dgemm(TotalMatrix, NT, TNT, ((MNStruct *)globalcontext)->pulse->nobs,
                 totalsize, ((MNStruct *)globalcontext)->pulse->nobs, totalsize,
                 'T', 'N');

    vector_dgemv(NT, Resvec, NTd, ((MNStruct *)globalcontext)->pulse->nobs,
                 totalsize, 'T');
  }
  if (savememory == 1) {
    for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
      for (int j = 0; j < totalsize; j++) {
        TotalMatrix[i + j * ((MNStruct *)globalcontext)->pulse->nobs] *=
            std::sqrt(Noise[i]);
      }
      Resvec[i] *= Noise[i];
    }

    vector_dgemm(TotalMatrix, TotalMatrix, TNT,
                 ((MNStruct *)globalcontext)->pulse->nobs, totalsize,
                 ((MNStruct *)globalcontext)->pulse->nobs, totalsize, 'T', 'N');

    for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
      for (int j = 0; j < totalsize; j++) {
        TotalMatrix[i + j * ((MNStruct *)globalcontext)->pulse->nobs] /=
            std::sqrt(Noise[i]);
      }
    }

    vector_dgemv(TotalMatrix, Resvec, NTd,
                 ((MNStruct *)globalcontext)->pulse->nobs, totalsize, 'T');
  }
  logtchk("Fnishing main algebra");

  if (savememory == 0) {
    delete[] NT;
  }
  for (int j = 0; j < totCoeff; j++) {
    TNT[TimetoMargin + j + (TimetoMargin + j) * totalsize] +=
        1.0 / powercoeff[j];
  }
  for (int j = 0; j < totalredshapecoeff; j++) {
    freqdet = freqdet - log(pow(10.0, -12) *
                            TNT[TimetoMargin + totCoeff + j +
                                (TimetoMargin + totCoeff + j) * totalsize]);
    TNT[TimetoMargin + totCoeff + j +
        (TimetoMargin + totCoeff + j) * totalsize] +=
        pow(10.0, -12) * TNT[TimetoMargin + totCoeff + j +
                             (TimetoMargin + totCoeff + j) * totalsize];
  }

  double freqlike = 0;
  double *WorkCoeff = new double[totalsize];
  double *WorkCoeff2 = new double[totalsize];
  double *TNT2 = new double[totalsize * totalsize];

  std::memcpy(TNT2, TNT, totalsize * totalsize * sizeof(double));
  std::memcpy(WorkCoeff, NTd, totalsize * sizeof(double));
  std::memcpy(WorkCoeff2, NTd, totalsize * sizeof(double));

  int globalinfo = 0;
  int info = 0;
  double jointdet = 0;
  vector_dpotrfInfo(TNT, totalsize, jointdet, info);
  if (info != 0) globalinfo = 1;

  info = 0;
  vector_dpotrsInfo(TNT, WorkCoeff, totalsize, info);

  if (info != 0) globalinfo = 1;
  info = 0;
  double det2 = 0;
  vector_TNqrsolve(TNT2, NTd, WorkCoeff2, totalsize, det2, info);

  double freqlike2 = 0;
  for (int j = 0; j < totalsize; j++) {
    freqlike += NTd[j] * WorkCoeff[j];
    freqlike2 += NTd[j] * WorkCoeff2[j];
  }

  double lnewChol = -0.5 * (tdet + jointdet + freqdet + timelike - freqlike) +
                    uniformpriorterm;
  double lnew =
      -0.5 * (tdet + det2 + freqdet + timelike - freqlike2) + uniformpriorterm;

  std::cout << "lnew " << lnew << " " << tdet << " " << jointdet << " "
            << freqdet << " " << timelike << " " << freqlike << " "
            << uniformpriorterm << std::endl;
  // printf("Double %g %i\n", lnew, globalinfo);
  if (fabs(lnew - lnewChol) > 0.05) {
    globalinfo = 1;
    //	lnew=-pow(10.0,20);
  }
  if (isnan(lnew) || isinf(lnew) || globalinfo != 0) {
    globalinfo = 1;
    lnew = -pow(10.0, 20);
  }

  if (badshape == 1) {
    lnew = -pow(10.0, 20);
  }

  delete[] DMVec;
  delete[] ScatVec;
  delete[] WorkCoeff;
  delete[] WorkCoeff2;
  for (int i = 0; i < ((MNStruct *)globalcontext)->EPolTerms; i++)
    delete[] EFAC[i];
  delete[] EFAC;
  delete[] EQUAD;
  delete[] SQUAD;
  delete[] powercoeff;
  delete[] NTd;
  delete[] freqs;
  delete[] Noise;
  delete[] Resvec;
  delete[] TNT;
  delete[] TNT2;

  if (((MNStruct *)globalcontext)->incNGJitter > 0) {
    delete[] ECORRPrior;
  }

  if (((MNStruct *)globalcontext)->incNGSJitter > 0) {
    delete[] SECORRPrior;
  }

  if (totalredshapecoeff > 0) {
    for (int j = 0; j < ((MNStruct *)globalcontext)->pulse->nobs; j++) {
      delete[] RedShapeletMatrix[j];
    }
    delete[] RedShapeletMatrix;
  }
  // printf("brace? %i \n", ((MNStruct *)globalcontext)->pulse->brace);
  if (((MNStruct *)globalcontext)->pulse->brace == -1) {
    lnew = -pow(10.0, 10);
  }
  ((MNStruct *)globalcontext)->PreviousInfo = globalinfo;
  ((MNStruct *)globalcontext)->PreviousJointDet = jointdet;
  ((MNStruct *)globalcontext)->PreviousFreqDet = freqdet;
  ((MNStruct *)globalcontext)->PreviousUniformPrior = uniformpriorterm;

  // printf("tdet %g, jointdet %g, freqdet %g, lnew %g, timelike %g, freqlike
  // %g\n", tdet, jointdet, freqdet, lnew, timelike, freqlike);

  // printf("CPUChisq: %g %g %g %g %g %g
  // \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
  logtchk("Exiting TempoNest Likelihood");

  return lnew;
}

void TemplateProfLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew,
                            void *context) {
  for (int p = 0; p < ndim; p++) {
    Cube[p] = (((MNStruct *)globalcontext)->PriorsArray[p + ndim] -
               ((MNStruct *)globalcontext)->PriorsArray[p]) *
                  Cube[p] +
              ((MNStruct *)globalcontext)->PriorsArray[p];
  }

  double *DerivedParams = new double[npars];

  double result = TemplateProfLike(ndim, Cube, npars, DerivedParams, context);

  delete[] DerivedParams;

  lnew = result;
}

double TemplateProfLike(int &ndim, double *Cube, int &npars,
                        double *DerivedParams, void *context) {
  // printf("In like \n");

  double *EFAC;
  double *EQUAD;
  double TemplateFlux = 0;
  double FakeRMS = 0;
  double TargetSN = 0;
  long double LDparams;
  int pcount;

  double uniformpriorterm = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Timing
  /// Model////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  LDparams = Cube[0] * (((MNStruct *)globalcontext)->LDpriors[0][1]) +
             (((MNStruct *)globalcontext)->LDpriors[0][0]);

  pcount = 0;
  long double phase =
      LDparams * ((MNStruct *)globalcontext)->ReferencePeriod / SECDAY;
  pcount++;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Profile
  /// Params//////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  long double **ProfileBats =
      new long double *[((MNStruct *)globalcontext)->pulse->nobs];
  long double *ModelBats =
      new long double[((MNStruct *)globalcontext)->pulse->nobs];
  TemplateFlux = 0;
  for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
    int Nbin = (int)((MNStruct *)globalcontext)->ProfileInfo[i][1];
    ProfileBats[i] = new long double[Nbin];
    for (int j = 0; j < Nbin; j++) {
      ProfileBats[i][j] = ((MNStruct *)globalcontext)->ProfileData[i][j][0];
    }

    ModelBats[i] = ((MNStruct *)globalcontext)->ProfileData[i][0][0] + phase;
  }

  int totalProfCoeff = 0;
  int *numcoeff = new int[((MNStruct *)globalcontext)->numProfComponents];
  // numcoeff[0] = ((MNStruct *)globalcontext)->numshapecoeff;
  // numcoeff[1] = 10;
  for (int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
    numcoeff[i] = ((MNStruct *)globalcontext)->numshapecoeff[i];
    totalProfCoeff += numcoeff[i];
  }

  double beta = Cube[pcount] * ((MNStruct *)globalcontext)->ReferencePeriod;

  pcount++;

  long double *CompSeps =
      new long double[((MNStruct *)globalcontext)->numProfComponents];
  CompSeps[0] = 0;
  for (int i = 1; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
    CompSeps[i] = Cube[pcount] * ((MNStruct *)globalcontext)->ReferencePeriod;
    pcount++;
    // printf("CompSep: %g \n", Cube[pcount-1]);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Profiles////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double lnew = 0;
  int badshape = 0;

  for (int t = 0; t < ((MNStruct *)globalcontext)->pulse->nobs; t++) {
    //	printf("In toa %i \n", t);
    int nTOA = t;

    double detN = 0;
    double logMargindet = 0;
    double Marginlike = 0;

    double profilelike = 0;

    long double FoldingPeriod =
        ((MNStruct *)globalcontext)->ProfileInfo[nTOA][0];
    long double FoldingPeriodDays = FoldingPeriod / SECDAY;
    int Nbins = (int)((MNStruct *)globalcontext)->ProfileInfo[nTOA][1];
    double Tobs = (double)((MNStruct *)globalcontext)->ProfileInfo[nTOA][2];
    long double ReferencePeriod = ((MNStruct *)globalcontext)->ReferencePeriod;

    double *Betatimes = new double[Nbins];
    double **Hermitepoly = new double *[Nbins];

    for (int i = 0; i < Nbins; i++) {
      Hermitepoly[i] = new double[totalProfCoeff];
      for (int j = 0; j < totalProfCoeff; j++) {
        Hermitepoly[i][j] = 0;
      }
    }

    long double binpos = ModelBats[nTOA];

    if (binpos < ProfileBats[nTOA][0]) binpos += FoldingPeriodDays;

    long double minpos = binpos - FoldingPeriodDays / 2;
    if (minpos < ProfileBats[nTOA][0]) minpos = ProfileBats[nTOA][0];
    long double maxpos = binpos + FoldingPeriodDays / 2;
    if (maxpos > ProfileBats[nTOA][Nbins - 1])
      maxpos = ProfileBats[nTOA][Nbins - 1];

    int cpos = 0;
    for (int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
      for (int j = 0; j < Nbins; j++) {
        long double timediff = 0;
        long double bintime = ProfileBats[t][j] + CompSeps[i] / SECDAY;
        if (bintime >= minpos && bintime <= maxpos) {
          timediff = bintime - binpos;
        } else if (bintime < minpos) {
          timediff = FoldingPeriodDays + bintime - binpos;
        } else if (bintime > maxpos) {
          timediff = bintime - FoldingPeriodDays - binpos;
        }

        timediff = timediff * SECDAY;

        Betatimes[j] = (timediff) / beta;
        TNothplMC(numcoeff[i], Betatimes[j], Hermitepoly[j], cpos);

        for (int k = 0; k < numcoeff[i]; k++) {
          // if(k==0)printf("HP %i %i %g %g \n", i, j,
          // Betatimes[j],Hermitepoly[j][cpos+k]*exp(-0.5*Betatimes[j]*Betatimes[j]));
          double Bconst =
              (1.0 / std::sqrt(beta)) /
              std::sqrt(
                  std::pow(2.0, k) *
                  std::sqrt(
                      M_PI));  //*((MNStruct *)globalcontext)->Factorials[k]);
          Hermitepoly[j][cpos + k] = Hermitepoly[j][cpos + k] * Bconst *
                                     exp(-0.5 * Betatimes[j] * Betatimes[j]);
        }
      }
      cpos += numcoeff[i];
    }

    ///////////////////////////////////////////Marginalise over arbitrary offset
    /// and absolute
    /// amplitude////////////////////////////////////////////////////////////

    double **M = new double *[Nbins];

    int Msize = totalProfCoeff + 1;
    for (int i = 0; i < Nbins; i++) {
      M[i] = new double[Msize];

      M[i][0] = 1;

      for (int j = 0; j < totalProfCoeff; j++) {
        M[i][j + 1] = Hermitepoly[i][j];
        // if(j==0)printf("%i %i %g \n", i, j, M[i][j+1]);
      }
    }

    double **MNM = new double *[Msize];
    for (int i = 0; i < Msize; i++) {
      MNM[i] = new double[Msize];
    }

    dgemm(M, M, MNM, Nbins, Msize, Nbins, Msize, 'T', 'N');

    for (int j = 0; j < Msize; j++) {
      MNM[j][j] += pow(10.0, -14);
    }

    double minFlux = pow(10.0, 100);
    for (int j = 0; j < Nbins; j++) {
      if ((double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] <
          minFlux) {
        minFlux = (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1];
      }
    }

    double *dNM = new double[Msize];
    double *TempdNM = new double[Msize];

    double *NDiffVec = new double[Nbins];
    TemplateFlux = 0;
    for (int j = 0; j < Nbins; j++) {
      TemplateFlux +=
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          minFlux;
      NDiffVec[j] =
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          minFlux;
    }
    dgemv(M, NDiffVec, dNM, Nbins, Msize, 'T');

    for (int i = 0; i < Msize; i++) {
      TempdNM[i] = dNM[i];
    }

    int info = 0;
    double Margindet = 0;
    dpotrfInfo(MNM, Msize, Margindet, info);
    dpotrs(MNM, TempdNM, Msize);

    //		printf("Margindet %g \n", Margindet);
    double *shapevec = new double[Nbins];
    dgemv(M, TempdNM, shapevec, Nbins, Msize, 'N');

    TargetSN = 10000.0;
    FakeRMS = TemplateFlux / TargetSN;
    FakeRMS = 1.0 / (FakeRMS * FakeRMS);
    double Chisq = 0;
    for (int j = 0; j < Nbins; j++) {
      double datadiff =
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          shapevec[j] - minFlux;
      Chisq += datadiff * datadiff * FakeRMS;
      //				printf("%i %.10g %.10g \n", j,
      //(double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] ,
      // shapevec[j]);
    }

    profilelike = -0.5 * (Chisq);
    lnew += profilelike;

    delete[] shapevec;
    delete[] NDiffVec;
    delete[] dNM;
    delete[] Betatimes;

    for (int j = 0; j < Nbins; j++) {
      delete[] Hermitepoly[j];
      delete[] M[j];
    }
    delete[] Hermitepoly;
    delete[] M;

    for (int j = 0; j < Msize; j++) {
      delete[] MNM[j];
    }
    delete[] MNM;
  }

  for (int j = 0; j < ((MNStruct *)globalcontext)->pulse->nobs; j++) {
    delete[] ProfileBats[j];
  }
  delete[] ProfileBats;
  delete[] ModelBats;
  delete[] numcoeff;
  delete[] CompSeps;

  // printf("End Like: %.10g %g %g\n", lnew, TemplateFlux, FakeRMS);

  return lnew;
}

void WriteMaxTemplateProf(std::string longname, int &ndim) {
  // printf("In like \n");

  double *Cube = new double[ndim];
  int number_of_lines = 0;

  std::ifstream checkfile;
  std::string checkname;
  if (((MNStruct *)globalcontext)->sampler == 0) {
    checkname = longname + "phys_live.points";
  }
  if (((MNStruct *)globalcontext)->sampler == 1) {
    checkname = longname + "_phys_live.txt";
  }
  printf("%i %s \n", ((MNStruct *)globalcontext)->sampler, checkname.c_str());
  checkfile.open(checkname.c_str());
  std::string line;
  while (getline(checkfile, line)) ++number_of_lines;

  checkfile.close();

  std::ifstream summaryfile;
  std::string fname;
  if (((MNStruct *)globalcontext)->sampler == 0) {
    fname = longname + "phys_live.points";
  }
  if (((MNStruct *)globalcontext)->sampler == 1) {
    fname = longname + "_phys_live.txt";
  }
  summaryfile.open(fname.c_str());

  printf("Getting ML \n");
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
        Cube[i] = paramlist[i];
      }
    }
  }
  summaryfile.close();

  double *EFAC;
  double *EQUAD;
  double TemplateFlux = 0;
  double FakeRMS = 0;
  double TargetSN = 0;
  long double LDparams;
  int pcount;

  double uniformpriorterm = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Timing
  /// Model////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  LDparams = Cube[0] * (((MNStruct *)globalcontext)->LDpriors[0][1]) +
             (((MNStruct *)globalcontext)->LDpriors[0][0]);

  pcount = 0;
  long double phase =
      LDparams * ((MNStruct *)globalcontext)->ReferencePeriod / SECDAY;
  pcount++;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Profile
  /// Params//////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  long double **ProfileBats =
      new long double *[((MNStruct *)globalcontext)->pulse->nobs];
  long double *ModelBats =
      new long double[((MNStruct *)globalcontext)->pulse->nobs];
  TemplateFlux = 0;
  for (int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++) {
    int Nbin = (int)((MNStruct *)globalcontext)->ProfileInfo[i][1];
    ProfileBats[i] = new long double[Nbin];
    for (int j = 0; j < Nbin; j++) {
      ProfileBats[i][j] = ((MNStruct *)globalcontext)->ProfileData[i][j][0];
    }

    ModelBats[i] = ((MNStruct *)globalcontext)->ProfileData[i][0][0] + phase;
  }

  int totalProfCoeff = 0;
  int *numcoeff = new int[((MNStruct *)globalcontext)->numProfComponents];

  for (int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
    numcoeff[i] = ((MNStruct *)globalcontext)->numshapecoeff[i];
    totalProfCoeff += numcoeff[i];
  }

  double beta = Cube[pcount] * ((MNStruct *)globalcontext)->ReferencePeriod;

  pcount++;

  long double *CompSeps =
      new long double[((MNStruct *)globalcontext)->numProfComponents];
  CompSeps[0] = 0;
  for (int i = 1; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
    CompSeps[i] = Cube[pcount] * ((MNStruct *)globalcontext)->ReferencePeriod;
    pcount++;
    // printf("CompSep: %g \n", Cube[pcount-1]);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////Profiles////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  double lnew = 0;
  int badshape = 0;

  for (int t = 0; t < ((MNStruct *)globalcontext)->pulse->nobs; t++) {
    //	printf("In toa %i \n", t);
    int nTOA = t;

    double detN = 0;
    double logMargindet = 0;
    double Marginlike = 0;

    double profilelike = 0;

    long double FoldingPeriod =
        ((MNStruct *)globalcontext)->ProfileInfo[nTOA][0];
    long double FoldingPeriodDays = FoldingPeriod / SECDAY;
    int Nbins = (int)((MNStruct *)globalcontext)->ProfileInfo[nTOA][1];
    double Tobs = (double)((MNStruct *)globalcontext)->ProfileInfo[nTOA][2];
    long double ReferencePeriod = ((MNStruct *)globalcontext)->ReferencePeriod;

    double *Betatimes = new double[Nbins];
    double **Hermitepoly = new double *[Nbins];

    for (int i = 0; i < Nbins; i++) {
      Hermitepoly[i] = new double[totalProfCoeff];
      for (int j = 0; j < totalProfCoeff; j++) {
        Hermitepoly[i][j] = 0;
      }
    }

    long double binpos = ModelBats[nTOA];

    if (binpos < ProfileBats[nTOA][0]) binpos += FoldingPeriodDays;

    long double minpos = binpos - FoldingPeriodDays / 2;
    if (minpos < ProfileBats[nTOA][0]) minpos = ProfileBats[nTOA][0];
    long double maxpos = binpos + FoldingPeriodDays / 2;
    if (maxpos > ProfileBats[nTOA][Nbins - 1])
      maxpos = ProfileBats[nTOA][Nbins - 1];

    int cpos = 0;
    for (int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
      printf("Sep %.10Lg \n", CompSeps[i]);
      for (int j = 0; j < Nbins; j++) {
        long double timediff = 0;
        long double bintime = ProfileBats[t][j] + CompSeps[i] / SECDAY;
        if (bintime >= minpos && bintime <= maxpos) {
          timediff = bintime - binpos;
        } else if (bintime < minpos) {
          timediff = FoldingPeriodDays + bintime - binpos;
        } else if (bintime > maxpos) {
          timediff = bintime - FoldingPeriodDays - binpos;
        }

        timediff = timediff * SECDAY;

        Betatimes[j] = (timediff) / beta;
        TNothplMC(numcoeff[i], Betatimes[j], Hermitepoly[j], cpos);

        for (int k = 0; k < numcoeff[i]; k++) {
          // if(k==0)printf("HP %i %i %g %g \n", i, j,
          // Betatimes[j],Hermitepoly[j][cpos+k]*exp(-0.5*Betatimes[j]*Betatimes[j]));
          double Bconst =
              (1.0 / std::sqrt(beta)) /
              std::sqrt(
                  pow(2.0, k) *
                  std::sqrt(
                      M_PI));  //*((MNStruct *)globalcontext)->Factorials[k]);
          Hermitepoly[j][cpos + k] = Hermitepoly[j][cpos + k] * Bconst *
                                     exp(-0.5 * Betatimes[j] * Betatimes[j]);
        }
      }
      cpos += numcoeff[i];
    }

    ///////////////////////////////////////////Marginalise over arbitrary offset
    /// and absolute
    /// amplitude////////////////////////////////////////////////////////////

    double **M = new double *[Nbins];

    int Msize = totalProfCoeff + 1;
    for (int i = 0; i < Nbins; i++) {
      M[i] = new double[Msize];

      M[i][0] = 1;

      for (int j = 0; j < totalProfCoeff; j++) {
        M[i][j + 1] = Hermitepoly[i][j];
        // if(j==0)printf("%i %i %g \n", i, j, M[i][j+1]);
      }
    }

    double **MNM = new double *[Msize];
    for (int i = 0; i < Msize; i++) {
      MNM[i] = new double[Msize];
    }

    dgemm(M, M, MNM, Nbins, Msize, Nbins, Msize, 'T', 'N');

    for (int j = 0; j < Msize; j++) {
      MNM[j][j] += pow(10.0, -14);
    }

    double minFlux = pow(10.0, 100);
    for (int j = 0; j < Nbins; j++) {
      if ((double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] <
          minFlux) {
        minFlux = (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1];
      }
    }

    double *dNM = new double[Msize];
    double *TempdNM = new double[Msize];

    double *NDiffVec = new double[Nbins];
    TemplateFlux = 0;
    for (int j = 0; j < Nbins; j++) {
      TemplateFlux +=
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          minFlux;
      NDiffVec[j] =
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          minFlux;
    }
    dgemv(M, NDiffVec, dNM, Nbins, Msize, 'T');

    for (int i = 0; i < Msize; i++) {
      TempdNM[i] = dNM[i];
    }

    int info = 0;
    double Margindet = 0;
    dpotrfInfo(MNM, Msize, Margindet, info);
    dpotrs(MNM, TempdNM, Msize);

    //		printf("Margindet %g \n", Margindet);
    double *shapevec = new double[Nbins];
    dgemv(M, TempdNM, shapevec, Nbins, Msize, 'N');

    TargetSN = 10000.0;
    FakeRMS = TemplateFlux / TargetSN;
    FakeRMS = 1.0 / (FakeRMS * FakeRMS);
    double Chisq = 0;
    for (int j = 0; j < Nbins; j++) {
      double datadiff =
          (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] -
          shapevec[j] - minFlux;
      Chisq += datadiff * datadiff * FakeRMS;
      printf("Bin %i Data %.10g Template %.10g Noise %.10g \n", j,
             (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1],
             shapevec[j], 1.0 / std::sqrt(FakeRMS));
    }

    cpos = 0;
    for (int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++) {
      for (int j = 0; j < numcoeff[i]; j++) {
        printf("P%i %.10g \n", i, TempdNM[cpos + j + 1] / TempdNM[1]);
      }
      cpos += numcoeff[i];
    }
    printf("B %.10g \n", Cube[1]);

    profilelike = -0.5 * (Chisq);
    lnew += profilelike;

    delete[] shapevec;
    delete[] NDiffVec;
    delete[] dNM;
    delete[] Betatimes;

    for (int j = 0; j < Nbins; j++) {
      delete[] Hermitepoly[j];
      delete[] M[j];
    }
    delete[] Hermitepoly;
    delete[] M;

    for (int j = 0; j < Msize; j++) {
      delete[] MNM[j];
    }
    delete[] MNM;
  }

  for (int j = 0; j < ((MNStruct *)globalcontext)->pulse->nobs; j++) {
    delete[] ProfileBats[j];
  }
  delete[] ProfileBats;
  delete[] ModelBats;
  delete[] numcoeff;
  delete[] CompSeps;

  // printf("End Like: %.10g %g %g\n", lnew, TemplateFlux, FakeRMS);
}

int llongturn_hms(long double turn, char *hms) {
  /* Converts double turn to string " hh:mm:ss.ssss" */

  int hh, mm, isec;
  long double sec;

  hh = (int)(turn * 24.);
  mm = (int)((turn * 24. - hh) * 60.);
  sec = ((turn * 24. - hh) * 60. - mm) * 60.;
  isec = (int)((sec * 10000. + 0.5) / 10000);
  if (isec == 60) {
    sec = 0.;
    mm = mm + 1;
    if (mm == 60) {
      mm = 0;
      hh = hh + 1;
      if (hh == 24) {
        hh = 0;
      }
    }
  }

  sprintf(hms, " %02d:%02d:%.12Lg", hh, mm, sec);
  return 0;
}

int llongturn_dms(long double turn, char *dms) {
  /* Converts double turn to string "sddd:mm:ss.sss" */

  int dd, mm, isec;
  long double trn, sec;
  char sign;

  sign = ' ';
  if (turn < 0.) {
    sign = '-';
    trn = -turn;
  } else {
    sign = '+';
    trn = turn;
  }
  dd = (int)(trn * 360.);
  mm = (int)((trn * 360. - dd) * 60.);
  sec = ((trn * 360. - dd) * 60. - mm) * 60.;
  isec = (int)((sec * 1000. + 0.5) / 1000);
  if (isec == 60) {
    sec = 0.;
    mm = mm + 1;
    if (mm == 60) {
      mm = 0;
      dd = dd + 1;
    }
  }
  sprintf(dms, "%c%02d:%02d:%.12Lg", sign, dd, mm, sec);
  return 0;
}
