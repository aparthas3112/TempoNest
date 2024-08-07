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
    double Tspan;
    int TimetoMargin;
    int totRedShapeCoeff;
    int totalsize;

    char* rootName;
    long double** LDpriors;
    double** Dpriors;
    double** DMatrix;
    double** FMatrix;
    double** GMatrix;
    double**
        staticGMatrix;  // staticGMatrix and staicTimedet for speeding up when no fitted white noise
    double staticTimeDet;
    double** UMatrix;  // UMatrix and SVec are for speeding up when only have 1EFAC and 1EQUAD
    double* SVec;
    double* maxLikeRes;
    int numberpulsars;
    int doLinear;
    int numFitJumps;
    int numFitTiming;
    int numFitEFAC;

    int numFitEQUAD;
    int incDMEQUAD;
    int* includeEQsys;
    int numFitRedCoeff;
    int numFitDMCoeff;

    int totCoeff;
    int numFitRedPL;
    int numFitDMPL;

    double* sampleFreq;
    int numdims;
    int incRED;
    int incDM;

    int Gsize;
    int Dsize;
    int** TempoFitNums;
    int* TempoJumpNums;
    int* sysFlags;
    double* TobsInfo;
    int systemcount;
    int TimeMargin;
    int JumpMargin;
    std::string* name;

    char* whiteflag;

    int whitemodel;

    int RedPriorType;
    int DMPriorType;

    int EQUADPriorType;
    int EFACPriorType;
    int useOriginalErrors;

    int FitLowFreqCutoff;
    int useNbitsAlgebra;

    int printResiduals;

    int storeFMatrices;
    double* StoredTMatrix;
    double* StoredDMVals;
    int* TimingGradientSigns;
    int usecosiprior;
    double* PreJumpVals;
    int doMax;

    double* DMatrixVec;
    double* PriorsArray;

    int debug;
    int rank;
} MNStruct;

void assigncontext(void* context);
void assignGPUcontext(void* context);
void assignGHScontext(void* context);

void callGHS(int NBurn, int NSamp, int GHSresume);
void callFourierDomainGHS(int NBurn, int NSamp, int GHSresume);

double iter_factorial(unsigned int n);
void store_factorial();
void fastephemeris_routines(pulsar* psr, int npsr);
void fastSubIntephemeris_routines(pulsar* psr, int npsr);
void fastformBatsAll(pulsar* psr, int npsr);
void fastformSubIntBatsAll(pulsar* psr, int npsr);
void outputProfile(int ndim);

void TNtextOutput(pulsar* psr, int npsr, int newpar, long double* Tempo2Fit, void* context,
                  int incRED, int ndims, std::vector<double> paramlist, double Evidence,
                  int MarginTime, int MarginJumps, int doLinear, std::string longname,
                  double** paramarray);
void getmaxlikeDM(pulsar* pulse, std::string longname, int ndim, void* context,
                  double** paramsarray);
void getProfileNoiseLevels(void* context);

double AllTOALike(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
double AllTOAJitterLike(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
double AllTOASim(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
double AllTOAStocProfLike(int& ndim, double* Cube, int& npars, double* DerivedParams,
                          void* context);
double AllTOAMaxLike(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
double AllTOAWriteMaxLike(std::string longname, int& ndim);
double AllTOAMarginStocProfLike(int& ndim, double* Cube, int& npars, double* DerivedParams,
                                void* context);
double TemplateProfLike(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
void WriteMaxTemplateProf(std::string longname, int& ndim);
void WriteSubIntStocProfLike(std::string longname, int& ndim);
void PreComputeShapelets(double** StoredShapelet, double** StoredJitter, double** StoredWidth,
                         double** InterpolatedMeanProfile, double** InterpolatedJitterProfile,
                         double** InterpolatedWidthProfile, long double finalInterpTime,
                         int numtointerpolate, double* MeanBeta, double& MaxShapeAmp);
void Tscrunch(void* globalcontext, double TemplateChanWidth);
void Tscrunch2(void* globalcontext, double TemplateChanWidth);
void getNumTempFreqs(int& NumFreqs, void* context, double TemplateChanWidth);
// void FFTProfileData(double **StoredFourierShapelets, double **StoredFourierJitterShapelets, int
// &NFBasis, int FindThreshold, double **FourierProfileData, void *context);

double ProfileDomainLike(int& ndim, double* Cube, int& npars, double* DerivedParams, void* context);
void ProfileDomainLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
void WriteProfileDomainLike(std::string longname, int& ndim);

double Template2DProfLike(int& ndim, double* Cube, int& npars, double* DerivedParams,
                          void* context);
void Template2DProfLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
void WriteTemplate2DProfLike(std::string longname, int& ndim);
void WriteDMSignal(std::string longname, int& ndim);
void TemplateProfLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
double FitPhase();
void GetInterpVectors(void* context);

double SubIntStocProfLike(int& ndim, double* Cube, int& npars, double* DerivedParams,
                          void* context);
void SubIntStocProfLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
// non linear timing model likelihood functions
// double  WhiteLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
// double NewLRedMarginLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void
// *context);
double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived, void* context);
// void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
// double LRedNumericalLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void
// *context);

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
void LRedGPULikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);
// non-Gaussian noise likelihoods
void subtractMLsolution(void* context);
void processPDF(void* context, std::string longname);
void TempoNestNGConvolvedLikeFunc(double* Cube, int& ndim, int& npars, double& lnew, void* context);
void TempoNestNGLikeFunc(double* Cube, int& ndim, int& npars, double& lnew, void* context);

// GPU non linear timing model likelihood functions
// void WhiteMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
// void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double NewLRedMarginGPULogLike(int& ndim, double* Cube, int& npars, double* DerivedParams,
                               void* globalcontext);
void NelderMeadOptimum(int nParameters);

// Functions to calculate the design matrices or 'G' marginalisation matrices
// void makeGDesign(pulsar *pulse, int &Gsize, int numtofit, double** staticGMatrix, double
// **oneDesign); void getDMatrix(pulsar *pulse, int TimeToFit, int JumptoFit, int numToMargin, int
// **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double
// **TNDM); void getMarginDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int
// **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double
// **TNDM, int linearFit);
void getCustomDMatrix(pulsar* pulse, int* MarginList, int** TempoFitNums, int* TempoJumpNums,
                      double** Dpriors, int incDM, int TimetoFit, int JumptoFit);
// void makeStaticGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** staticGMatrix, double
// &tdet); void makeStaticDiagGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** UMatrix,
// double *SVec);
void getCustomDMatrixLike(void* context, double** TNDM);
void getCustomDVectorLike(void* context, double* TNDM, int nobs, int TimeToMargin, int TotalSize);

void getNGJitterMatrix(pulsar* pulse, double** JitterMatrix, int& NumEpochs);
void getNGSJitterMatrix(pulsar* pulse, double** JitterMatrix, int& NumEpochs);

void getNGJitterMatrixEpochs(pulsar* pulse, int& NumEpochs);
void getNGSJitterMatrixEpochs(pulsar* pulse, int& NumEpochs);

void StoreTMatrix(double* TMatrix, void* context);
void getArraySizeInfo(void* context);
void getPhysDVector(void* context, double** TNDM, int Nobs, int* TimingGradientSigns);
void UpdatePhysDVector(void* context, double** TNDM, int Nobs);
void OutputMLFiles(int nParameters, double* pdParameterEstimates, double MLike, int startDim);

void readsummary(pulsar* psr, std::string longname, int ndim, void* context, long double* Tempo2Fit,
                 int incRED, int ndims, int MarginTime, int MarginJumps, int doLinear);

void setupMNparams(char* ConfigFileName, int& sampler, int& IS, int& modal, int& ceff, int& nlive,
                   double& efr, int& sample, int& updInt, int& nClsPar, int& Nchords, int& NBurn,
                   int& NSamp, int& GHSresume);
void setupparams(char* ConfigFileName, int& useGPUS, char* root, int& numTempo2its,
                 int& doLinearFit, int& doMax, int& incEFAC, int& incEQUAD, int& incRED, int& incDM,
                 int& doTimeMargin, int& doJumpMargin, double& FitSig, int& customPriors,
                 double* EFACPrior, double* EQUADPrior, double* AlphaPrior, double* AmpPrior,
                 double* DMAlphaPrior, double* DMAmpPrior, double& numRedCoeff, double& numDMCoeff,
                 double& numScatCoeff, int& numRedPL, int& numDMPL, double& FourierSig,
                 char* whiteflag, int& whitemodel, int& RedPriorType, int& DMPriorType,
                 int& EQUADPriorType, int& EFACPriorType, int& useOriginalErrors,
                 int& FitLowFreqCutoff, int& useNbitsAlgebra, int& StoreTMatrix, int& debug);

void setTNPriors(char* ConfigFileName, double** Dpriors, long double** TempoPriors, int TPsize,
                 int DPsize);
// void setFrequencies(char *ConfigFileName, double *SampleFreq, int numRedfreqs, int numDMfreqs,
// int numRedLogFreqs, int numDMLogFreqs, double RedLowFreq, double DMLowFreq, double RedMidFreq,
// double DMMidFreq);
void setFrequencies(char* ConfigFileName, double* SampleFreq, int numRedfreqs, int numDMfreqs,
                    int numScatfreqs, int numRedLogFreqs, int numDMLogFreqs, int numScatLogFreqs,
                    double RedLowFreq, double DMLowFreq, double ScatLowFreq, double RedMidFreq,
                    double DMMidFreq, double ScatMidFreq);
void GetGroupsToFit(char* ConfigFileName, int incGroupNoise, int** FitForGroup, int incBandNoise,
                    int** FitForBand);
void setShapePriors(char* ConfigFileName, double** ShapePriors, double** BetaPrior, int numcoeff,
                    int numcomps);
void GetProfileFitInfo(char* ConfigFileName, int numProfComponents, int* numGPTAshapecoeff,
                       int* numProfileFitCoeff, int* numEvoCoeff, int* numFitEvoCoeff,
                       int* numGPTAstocshapecoeff, double* ProfCompSeps, double& TemplateChanWidth,
                       int* TimeCorrShapeCoeff, int incExtraComp, double** FitForExtraComp,
                       int* FitCompWidths, int* FitCompPos);

void TNprior(double cube[], double theta[], int nDims, void* context);
