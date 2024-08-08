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

#include <vector>
#include "tempo2.h"

typedef struct {
    pulsar* pulse;
    int numberpulsars;
    double Tspan;
    int TimetoMargin;
    int totalsize;

    char* rootName;
    double** Dpriors;

    int numFitJumps;
    int numFitTiming;
    int numFitEFAC;

    int numFitEQUAD;
    int* includeEQsys;
    int numFitRedCoeff;
    int numFitDMCoeff;

    int totCoeff;

    double* sampleFreq;
    int numdims;
    int incRED;
    int incDM;

    int** TempoFitNums;
    int* TempoJumpNums;
    int* sysFlags;
    int systemcount;
    int TimeMargin;
    int JumpMargin;

    char* whiteflag;

    int RedPriorType;
    int DMPriorType;

    int EQUADPriorType;
    int EFACPriorType;
    int useOriginalErrors;
    int printResiduals;

    int storeFMatrices;
    double* StoredTMatrix;
    int usecosiprior;
    double* PreJumpVals;

    double* DMatrixVec;
    double* PriorsArray;

    int debug;
    int rank;
} MNStruct;

void assigncontext(void* context);
void fastephemeris_routines(pulsar* psr, int npsr);
void fastSubIntephemeris_routines(pulsar* psr, int npsr);
void fastformBatsAll(pulsar* psr, int npsr);
void fastformSubIntBatsAll(pulsar* psr, int npsr);

void TNtextOutput(pulsar* psr, int npsr, int newpar, void* context, int incRED, int ndims,
                  std::vector<double> paramlist, double Evidence, std::string longname,
                  double** paramarray);

// non linear timing model likelihood functions
// double  WhiteLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
// double NewLRedMarginLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void
// *context);
double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived, void* context);
// void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
// double LRedNumericalLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void
// *context);

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);

void StoreTMatrix(double* TMatrix, void* context);
void getArraySizeInfo(void* context);
void OutputMLFiles(int nParameters, double* pdParameterEstimates, double MLike, int startDim);

void readsummary(pulsar* psr, std::string longname, int ndim, void* context, int incRED, int ndims);

void setupMNparams(char* ConfigFileName, int& sampler, int& IS, int& modal, int& ceff, int& nlive,
                   double& efr, int& sample, int& updInt, int& nClsPar, int& Nchords, int& NBurn,
                   int& NSamp, int& GHSresume);
void setupparams(char* ConfigFileName, char* root, int& numTempo2its, int& incEFAC, int& incEQUAD,
                 int& incRED, int& incDM, double* EFACPrior, double* EQUADPrior, double* AlphaPrior,
                 double* AmpPrior, double* DMAlphaPrior, double* DMAmpPrior, double& numRedCoeff,
                 double& numDMCoeff, double& FourierSig, char* whiteflag, int& RedPriorType,
                 int& DMPriorType, int& EQUADPriorType, int& EFACPriorType, int& useOriginalErrors,
                 int& StoreTMatrix, int& debug);

void setTNPriors(char* ConfigFileName, double** Dpriors, long double** TempoPriors, int TPsize,
                 int DPsize);
// void setFrequencies(char *ConfigFileName, double *SampleFreq, int numRedfreqs, int numDMfreqs,
// int numRedLogFreqs, int numDMLogFreqs, double RedLowFreq, double DMLowFreq, double RedMidFreq,
// double DMMidFreq);
void setFrequencies(char* ConfigFileName, double* SampleFreq, int numRedfreqs, int numDMfreqs,
                    int numRedLogFreqs, int numDMLogFreqs, int numScatLogFreqs, double RedLowFreq,
                    double DMLowFreq, double ScatLowFreq, double RedMidFreq, double DMMidFreq,
                    double ScatMidFreq);
