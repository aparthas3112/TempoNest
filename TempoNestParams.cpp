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
*    then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model and MultiNest Papers here.
*/


#include <string.h>
#include <stdio.h>
#include <math.h>
#include "configfile.h"

void setupparams(char *ConfigFileName, int &useGPUS,
		char *root,
		int &numTempo2its,
		int &doLinearFit, 
		int &doMax,
		int &incEFAC,
		int &EPolyTerms,
		int &incEQUAD,
		int &incRED,
		int &incDM,
		int &doTimeMargin,
		int &doJumpMargin,
		double &FitSig,
		int &customPriors,
		double *EFACPrior,
		double *EPolyPriors,
		double *EQUADPrior,
		double *AlphaPrior,
		double *AmpPrior,
		double *DMAlphaPrior,
		double *DMAmpPrior,
		double &numRedCoeff,
		double &numDMCoeff,
		int &numRedPL,
		int &numDMPL,
		double *RedCoeffPrior,
		double *DMCoeffPrior,
		int &FloatingDM,
		double *DMFreqPrior,
		int &yearlyDM,
		int &incsinusoid,
		int &FloatingRed,
		double *RedFreqPrior,
		double &FourierSig,
		int &incStep,
		double *StepAmpPrior,
		char *whiteflag,
		int &whitemodel,
		int &varyRedCoeff,
		int &varyDMCoeff,
		int &incGWB,
		double *GWBAmpPrior,
		int &RedPriorType,
		int &DMPriorType,
		int &EQUADPriorType,
		int &EFACPriorType,
		int &useOriginalErrors,
		int &incShannonJitter,		
		int &incDMEvent,
		double *DMEventStartPrior,
		double *DMEventLengthPrior,
		int &incDMShapeEvent,
		int &numDMShapeCoeff,
		double *DMShapeCoeffPrior,
		int &incRedShapeEvent,
		int &numRedShapeCoeff,
		int &MarginRedShapeCoeff,
		double *RedShapeCoeffPrior,
		int &incDMScatterShapeEvent,
		int &numDMScatterShapeCoeff,
		double *DMScatterShapeCoeffPrior,
		int &incBandNoise,
		int &numBandNoiseCoeff,
		double *BandNoiseAmpPrior,
		double *BandNoiseAlphaPrior,
		int &incNGJitter,
		int &incGlitch,
		int &incGlitchTerms,
		double &GlitchFitSig,
		int &incBreakingIndex,
		int &FitLowFreqCutoff,
		int &uselongdouble,
		int &incGroupNoise,
		int &numGroupCoeff,
		double *GroupNoiseAmpPrior,
		double *GroupNoiseAlphaPrior,
                int &FitSolarWind,
                int &FitWhiteSolarWind,
                double *SolarWindPrior,
                double *WhiteSolarWindPrior,
		int &GPTA,
		char *GroupNoiseFlag,
		int &FixProfile,
		int &FitTemplate,
		int &InterpolateProfile,
		double &InterpolatedTime,
		int &StoreTMatrix,
		int &incHighFreqStoc,
		double *HighFreqStocPrior,
		int &incProfileEvo,
		double &EvoRefFreq,
		double *ProfileEvoPrior,
		int &FitEvoExponent,
		int &incWideBandNoise,
		int &incProfileFit,
		double *ProfileFitPrior,
		int &FitLinearProfileWidth,
		double *LinearProfileWidthPrior,
		int &incDMEQUAD,
		double *DMEQUADPrior,
		double &offPulseLevel,
		char *ProfFile,
		int &numProfComponents,
		int &incWidthJitter,
		double *WidthJitterPrior,
		int &JitterProfComp,
		int &incProfileEnergyEvo,
		double *ProfileEnergyEvoPrior,
		int &debug,
		int &ProfileBaselineTerms,
                int &incProfileNoise,
                int &ProfileNoiseCoeff,
                double *ProfileNoiseAmpPrior,
                double *ProfileNoiseSpecPrior,
		int &SubIntToFit,
		int &ChannelToFit,
		int &NProfileEvoPoly,
		int &usecosiprior,
		int &incWidthEvoTime,
		int &incExtraProfComp,
		int &removeBaseline,
		int &incPrecession, 
		int &incTimeCorrProfileNoise,
		double &phasePriorExpansion,
		int &ProfileNoiseMethod,
		int &FitPrecAmps,
		int &NProfileTimePoly,
		int &incProfileScatter,
		int &ScatterPBF){

	//General parameters:
	//Use GPUs 0=No, 1=Yes
		
	debug = 0;
	useGPUS=0;
	uselongdouble=0;

	GPTA = 0;
	numProfComponents = 0;
	FixProfile = 0;	
	FitTemplate = 0;
	removeBaseline = 0;
	ProfileBaselineTerms = 1;

	InterpolateProfile = 0;
	InterpolatedTime = 1; //in nanoseconds
	phasePriorExpansion = 1;
	ProfileNoiseMethod = 0;


	StoreTMatrix = 1; // Recompute TMatrices when computing new bats - default is just precompute and use those



	SubIntToFit = -1;
	ChannelToFit = -1;

	//Root of the results files,relative to the directory in which TempoNest is run. This will be followed by the pulsar name, and then the individual output file extensions.
	strcpy( root, "results/Example1-");
	strcpy( ProfFile, "Profile.dat");
	//numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.  Should only be set to 0 if all the priors are set in setTNPriors
	numTempo2its=1;

	//doLinearFit:  Switches between the full non linear timing model (doLinearFit=0) and the linear approximation for the timing model based on the initial Tempo2 Fit (= 1).

	doLinearFit=0;

	//doMax: Find maximum likelihood values for non linear timing model for the stochastic model chosen.
	//Will find maximum in the full un marginalised problem in order to find the best values to then marginalise over.
	//Start point of non linear search will be performed at this maximum if chosen unless set otherwise in custom priors.
	//Central point of prior for non linear search will be set at this value unless set otherwise in custom priors.
	doMax=0;

	//ModelChoices

	//White noise

	useOriginalErrors = 0;  // Use tempo2 errors before modification by TNEF/TNEQ
	incEFAC=0; //include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
	EPolyTerms = 1; //Number of terms to include in EFAC polynomial (A*TOAerr + B*TOAERR^2 etc)
	incEQUAD=0; //include EQUAD: 0 = no, 1 = yes
	incNGJitter = 0; //Include NG Jitter for systems flagged in par file
	incDMEQUAD = 0;

	JitterProfComp = 0;

	FitSolarWind = 0; // Basically for for ne_sw
	FitWhiteSolarWind = 0; //Fit for a white component proportional to tdis2 with ne_sw=1


	strcpy( whiteflag, "-sys");
	whitemodel=0;

	incShannonJitter=0;

	incRED=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)
	incDM=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)

	FitLowFreqCutoff = 0; //Include f_low as a free parameter

	doTimeMargin=0 ; //0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
	doJumpMargin=0; //0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.
	incBreakingIndex=0; //Use breaking index to constrain F1/F2 etc


	incWideBandNoise = 0;
	EvoRefFreq = 1400.0;
	incProfileEvo = 0;
	NProfileEvoPoly = 1;
	FitEvoExponent = 0;
	ProfileEvoPrior[0] = -1;
	ProfileEvoPrior[1] =  1;

	incProfileEnergyEvo = 0;
	ProfileEnergyEvoPrior[0] = -1;
	ProfileEnergyEvoPrior[1] =  1;

	NProfileTimePoly = 0;

	incProfileFit = 0;
	
	ProfileFitPrior[0] = -1;
	ProfileFitPrior[1] =  1;

	FitLinearProfileWidth = 0;
	LinearProfileWidthPrior[0] = -0.01;
	LinearProfileWidthPrior[1] =  0.01;

	incWidthEvoTime = 0;

	incWidthJitter = 0;
	WidthJitterPrior[0] =  -10;
	WidthJitterPrior[1] =  -1;

	incProfileNoise = 0;
	ProfileNoiseCoeff = 10;
	ProfileNoiseAmpPrior[0] = -10;
	ProfileNoiseAmpPrior[1] = 0;
	ProfileNoiseSpecPrior[0] = -7;
	ProfileNoiseSpecPrior[1] = 7;
	
	incExtraProfComp = 0;
	incPrecession = 0;
	FitPrecAmps = -1;

	incTimeCorrProfileNoise = 0;

	incProfileScatter = 0;
	ScatterPBF = 1;

    //Priors

    //Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
    customPriors=0; 


	//FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
	//For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
	// For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
	FitSig=5;
	
	//Remaining priors for the stochastic parameters.  

	RedPriorType = 0; // 0 = Log, 1 = Uniform
	DMPriorType = 0;   // 0 = Log, 1 = Uniform
	EQUADPriorType = 0;   // 0 = Log, 1 = Uniform
	EFACPriorType = 0;   // 0 = Log, 1 = Uniform
	usecosiprior = 0; // 0 = uniform in sini, 1 = uniform in cosi


	EFACPrior[0]=0.1;
	EFACPrior[1]=10;

	EPolyPriors[0]=-20;
	EPolyPriors[1]=20;   //Prior on terms in EPoly > linear
	
	
	EQUADPrior[0]=-10;
	EQUADPrior[1]=-5;

	DMEQUADPrior[0] = -10;
	DMEQUADPrior[1] = -1;


	SolarWindPrior[0] = 0;
	SolarWindPrior[1] = 10;
	WhiteSolarWindPrior[0] = -3;
	WhiteSolarWindPrior[1] = 3;
	
	numRedPL=0;
	numDMPL=0;

	numRedCoeff=10;
	numDMCoeff=10;

	varyRedCoeff=0;
	varyDMCoeff=0;
	
	AlphaPrior[0]=1.1;
	AlphaPrior[1]=6.1;
	
	AmpPrior[0]=-20;
	AmpPrior[1]=-10;
	
	DMAlphaPrior[0]=1.1;
	DMAlphaPrior[1]=6.1;
	
	DMAmpPrior[0]=-18;
	DMAmpPrior[1]=-8;

	
	RedCoeffPrior[0]=-10;
	RedCoeffPrior[1]=0;
	
	DMCoeffPrior[0]=-10;
	DMCoeffPrior[1]=0;
	
	FourierSig = 5;


	
	FloatingDM = 0;
	DMFreqPrior[0]=1;
	DMFreqPrior[1]=100;
	
	FloatingRed = 0;
	RedFreqPrior[0]=1;
	RedFreqPrior[1]=100;

	yearlyDM=0;
	incsinusoid=0;
	
	incStep = 0;
	StepAmpPrior[0] = -1;
	StepAmpPrior[1] = 1;

	incGWB=0;
	GWBAmpPrior[0] = -20;
	GWBAmpPrior[1] = -12;


	incDMEvent = 0;

	incDMShapeEvent = 0;
	numDMShapeCoeff = 0;
	DMShapeCoeffPrior[0] = -0.01;
	DMShapeCoeffPrior[1] = 0.01;

	incRedShapeEvent = 0;
	numRedShapeCoeff = 0;
	MarginRedShapeCoeff = 0;
	RedShapeCoeffPrior[0] = -0.01;
	RedShapeCoeffPrior[1] = 0.01;

	incDMScatterShapeEvent = 0;
	numDMScatterShapeCoeff = 0;
	DMScatterShapeCoeffPrior[0] = -0.01;
	DMScatterShapeCoeffPrior[1] = 0.01;


	incBandNoise = 0;
	numBandNoiseCoeff = 10;
	BandNoiseAmpPrior[0] = -20;
	BandNoiseAmpPrior[1] = -8;
	BandNoiseAlphaPrior[0] = 0;
	BandNoiseAlphaPrior[1] = 7;


	strcpy( GroupNoiseFlag, "-sys");
	incGroupNoise = 0;
	numGroupCoeff = 10;
	GroupNoiseAmpPrior[0] = -20;
	GroupNoiseAmpPrior[1] = -8;
	GroupNoiseAlphaPrior[0] = 0;
	GroupNoiseAlphaPrior[1] = 7;



	//GPTA Params
	
	incHighFreqStoc = 0;
	HighFreqStocPrior[0] = -10;
	HighFreqStocPrior[1] = 1;
	
	offPulseLevel = 0.001;



    // Use a configfile, if we can, to overwrite the defaults set in this file.
    try {
        string strBuf;
	
        strBuf = ConfigFileName;//string("defaultparameters.conf");
	//string NewstrBuf = ConfigFileName;

        ConfigFile parameters(strBuf);

        /* We can check whether a value is not set in the file by doing
         * if(! parameters.readInto(variable, "name", default)) {
         *   printf("WARNING");
         * }
         *
         * At the moment I was too lazy to print warning messages, and the
         * default value from this file is used in that case.
         *
         * Note: the timing model parameters are not done implemented yet
         */
	parameters.readInto(useGPUS, "useGPUS", useGPUS);
	parameters.readInto(uselongdouble, "uselongdouble", uselongdouble);
	parameters.readInto(GPTA, "GPTA", GPTA);
	parameters.readInto(FixProfile, "FixProfile", FixProfile);
	parameters.readInto(FitTemplate, "FitTemplate", FitTemplate);
	parameters.readInto(InterpolateProfile, "InterpolateProfile", InterpolateProfile);
	parameters.readInto(InterpolatedTime, "InterpolatedTime", InterpolatedTime);
        parameters.readInto(StoreTMatrix, "StoreTMatrix",StoreTMatrix );



        parameters.readInto(strBuf, "root", string("results/Example1"));
        strcpy(root, strBuf.data());
        parameters.readInto(numTempo2its, "numTempo2its", numTempo2its);
        parameters.readInto(doLinearFit, "doLinearFit", doLinearFit);
        parameters.readInto(doMax, "doMax", doMax);
        parameters.readInto(incEFAC, "incEFAC", incEFAC);
        parameters.readInto(EPolyTerms, "EPolyTerms", EPolyTerms);
        parameters.readInto(incEQUAD, "incEQUAD", incEQUAD);


	parameters.readInto(FitSolarWind, "FitSolarWind", FitSolarWind);
        parameters.readInto(FitWhiteSolarWind, "FitWhiteSolarWind", FitWhiteSolarWind);
        parameters.readInto(SolarWindPrior[0], "SolarWindPrior[0]", SolarWindPrior[0]);
        parameters.readInto(SolarWindPrior[1], "SolarWindPrior[1]", SolarWindPrior[1]);
        parameters.readInto(WhiteSolarWindPrior[0], "WhiteSolarWindPrior[0]", WhiteSolarWindPrior[0]);
        parameters.readInto(WhiteSolarWindPrior[1], "WhiteSolarWindPrior[1]", WhiteSolarWindPrior[1]);




        parameters.readInto(incRED, "incRED", incRED);
		parameters.readInto(incDM, "incDM", incDM);
        parameters.readInto(doTimeMargin, "doTimeMargin", doTimeMargin);
        parameters.readInto(doJumpMargin, "doJumpMargin", doJumpMargin);
        parameters.readInto(customPriors, "customPriors", customPriors);
        parameters.readInto(FitSig, "FitSig", FitSig);
        parameters.readInto(EFACPrior[0], "EFACPrior[0]", EFACPrior[0]);
        parameters.readInto(EFACPrior[1], "EFACPrior[1]", EFACPrior[1]);
	parameters.readInto(EPolyPriors[0], "EPolyPriors[0]", EPolyPriors[0]);
	parameters.readInto(EPolyPriors[1], "EPolyPriors[1]", EPolyPriors[1]);
        parameters.readInto(EQUADPrior[0], "EQUADPrior[0]", EQUADPrior[0]);
        parameters.readInto(EQUADPrior[1], "EQUADPrior[1]", EQUADPrior[1]);
        parameters.readInto(AlphaPrior[0], "AlphaPrior[0]", AlphaPrior[0]);
        parameters.readInto(AlphaPrior[1], "AlphaPrior[1]", AlphaPrior[1]);
        parameters.readInto(AmpPrior[0], "AmpPrior[0]", AmpPrior[0]);
        parameters.readInto(AmpPrior[1], "AmpPrior[1]", AmpPrior[1]);
        parameters.readInto(numRedCoeff, "numRedCoeff", numRedCoeff);
		parameters.readInto(numDMCoeff, "numDMCoeff", numDMCoeff);
        parameters.readInto(numRedPL, "numRedPL", numRedPL);
		parameters.readInto(numDMPL, "numDMPL", numDMPL);
        parameters.readInto(RedCoeffPrior[0], "RedCoeffPrior[0]", RedCoeffPrior[0]);
        parameters.readInto(RedCoeffPrior[1], "RedCoeffPrior[1]", RedCoeffPrior[1]);
        parameters.readInto(DMCoeffPrior[0], "DMCoeffPrior[0]", DMCoeffPrior[0]);
        parameters.readInto(DMCoeffPrior[1], "DMCoeffPrior[1]", DMCoeffPrior[1]);
        parameters.readInto(DMAlphaPrior[0], "DMAlphaPrior[0]", DMAlphaPrior[0]);
        parameters.readInto(DMAlphaPrior[1], "DMAlphaPrior[1]", DMAlphaPrior[1]);
        parameters.readInto(DMAmpPrior[0], "DMAmpPrior[0]", DMAmpPrior[0]);
        parameters.readInto(DMAmpPrior[1], "DMAmpPrior[1]", DMAmpPrior[1]);
        
        parameters.readInto(FloatingDM, "FloatingDM", FloatingDM);
	parameters.readInto(yearlyDM, "yearlyDM", yearlyDM);
	parameters.readInto(incsinusoid, "incsinusoid", incsinusoid);
        parameters.readInto(DMFreqPrior[0], "DMFreqPrior[0]", DMFreqPrior[0]);
        parameters.readInto(DMFreqPrior[1], "DMFreqPrior[1]", DMFreqPrior[1]);
        parameters.readInto(FloatingRed, "FloatingRed", FloatingRed);
        parameters.readInto(RedFreqPrior[0], "RedFreqPrior[0]", RedFreqPrior[0]);
        parameters.readInto(RedFreqPrior[1], "RedFreqPrior[1]", RedFreqPrior[1]);

        parameters.readInto(incStep, "incStep", incStep);
        parameters.readInto(StepAmpPrior[0], "StepAmpPrior[0]", StepAmpPrior[0]);
        parameters.readInto(StepAmpPrior[1], "StepAmpPrior[1]", StepAmpPrior[1]);

	parameters.readInto(strBuf, "whiteflag", string("-sys"));
        strcpy(whiteflag, strBuf.data());
	parameters.readInto(whitemodel, "whitemodel", whitemodel);

        parameters.readInto(FourierSig, "FourierSig", FourierSig);

        parameters.readInto(varyRedCoeff, "varyRedCoeff", varyRedCoeff);
        parameters.readInto(varyDMCoeff, "varyDMCoeff", varyDMCoeff);
	parameters.readInto(incGWB, "incGWB", incGWB);
	parameters.readInto(GWBAmpPrior[0], "GWBAmpPrior[0]", GWBAmpPrior[0]);
	parameters.readInto(GWBAmpPrior[1], "GWBAmpPrior[1]", GWBAmpPrior[1]);
	parameters.readInto(RedPriorType, "RedPriorType", RedPriorType);
	parameters.readInto(DMPriorType, "DMPriorType", DMPriorType);
	parameters.readInto(EFACPriorType, "EFACPriorType", EFACPriorType);
	parameters.readInto(EQUADPriorType, "EQUADPriorType", EQUADPriorType);
	parameters.readInto(useOriginalErrors, "useOriginalErrors", useOriginalErrors);
	parameters.readInto(incShannonJitter, "incShannonJitter", incShannonJitter);
	parameters.readInto(incNGJitter, "incNGJitter", incNGJitter);
	parameters.readInto(incGlitch, "incGlitch", incGlitch);
        parameters.readInto(incGlitchTerms, "incGlitchTerms", incGlitchTerms);        
	parameters.readInto(GlitchFitSig, "GlitchFitSig", GlitchFitSig);	
	parameters.readInto(incBreakingIndex, "incBreakingIndex", incBreakingIndex);
	parameters.readInto(FitLowFreqCutoff, "FitLowFreqCutoff", FitLowFreqCutoff);

	parameters.readInto(incDMShapeEvent, "incDMShapeEvent", incDMShapeEvent);
	parameters.readInto(numDMShapeCoeff, "numDMShapeCoeff", numDMShapeCoeff);

	int n; char buffer [50];
	for (int i=0; i<incDMShapeEvent; i++) {
	  n=sprintf (buffer, "DMEventStartPrior[%i]", 2*i);
	  parameters.readInto(DMEventStartPrior[2*i], buffer, DMEventStartPrior[2*i]);
	  n=sprintf (buffer, "DMEventStartPrior[%i]", 2*i+1);
	  parameters.readInto(DMEventStartPrior[2*i+1], buffer, DMEventStartPrior[2*i+1]);

	  n=sprintf (buffer, "DMEventLengthPrior[%i]", 2*i);
	  parameters.readInto(DMEventLengthPrior[2*i], buffer, DMEventLengthPrior[2*i]);
          n=sprintf (buffer, "DMEventLengthPrior[%i]", 2*i+1);
          parameters.readInto(DMEventLengthPrior[2*i+1], buffer, DMEventLengthPrior[2*i+1]);
	}
	parameters.readInto(DMShapeCoeffPrior[0], "DMShapeCoeffPrior[0]", DMShapeCoeffPrior[0]);
	parameters.readInto(DMShapeCoeffPrior[1], "DMShapeCoeffPrior[1]", DMShapeCoeffPrior[1]);

	parameters.readInto(incRedShapeEvent, "incRedShapeEvent", incRedShapeEvent);
	parameters.readInto(numRedShapeCoeff, "numRedShapeCoeff", numRedShapeCoeff);
	parameters.readInto(MarginRedShapeCoeff, "MarginRedShapeCoeff", MarginRedShapeCoeff);
	parameters.readInto(RedShapeCoeffPrior[0], "RedShapeCoeffPrior[0]", RedShapeCoeffPrior[0]);
	parameters.readInto(RedShapeCoeffPrior[1], "RedShapeCoeffPrior[1]", RedShapeCoeffPrior[1]);

	parameters.readInto(incDMScatterShapeEvent, "incDMScatterShapeEvent", incDMScatterShapeEvent);
	parameters.readInto(numDMScatterShapeCoeff, "numDMScatterShapeCoeff", numDMScatterShapeCoeff);
	parameters.readInto(DMScatterShapeCoeffPrior[0], "DMScatterShapeCoeffPrior[0]", DMScatterShapeCoeffPrior[0]);
	parameters.readInto(DMScatterShapeCoeffPrior[1], "DMScatterShapeCoeffPrior[1]", DMScatterShapeCoeffPrior[1]);


	parameters.readInto(incBandNoise, "incBandNoise", incBandNoise);
	parameters.readInto(numBandNoiseCoeff, "numBandNoiseCoeff", numBandNoiseCoeff);
	parameters.readInto(BandNoiseAmpPrior[0], "BandNoiseAmpPrior[0]", BandNoiseAmpPrior[0]);
	parameters.readInto(BandNoiseAmpPrior[1], "BandNoiseAmpPrior[1]", BandNoiseAmpPrior[1]);
	parameters.readInto(BandNoiseAlphaPrior[0], "BandNoiseAlphaPrior[0]", BandNoiseAlphaPrior[0]);
	parameters.readInto(BandNoiseAlphaPrior[1], "BandNoiseAlphaPrior[1]", BandNoiseAlphaPrior[1]);


	parameters.readInto(strBuf, "GroupNoiseFlag", string("-sys"));
        strcpy(GroupNoiseFlag, strBuf.data());
	parameters.readInto(incGroupNoise, "incGroupNoise", incGroupNoise);
	parameters.readInto(numGroupCoeff, "numGroupCoeff", numGroupCoeff);
	parameters.readInto(GroupNoiseAmpPrior[0], "GroupNoiseAmpPrior[0]", GroupNoiseAmpPrior[0]);
	parameters.readInto(GroupNoiseAmpPrior[1], "GroupNoiseAmpPrior[1]", GroupNoiseAmpPrior[1]);
	parameters.readInto(GroupNoiseAlphaPrior[0], "GroupNoiseAlphaPrior[0]", GroupNoiseAlphaPrior[0]);
	parameters.readInto(GroupNoiseAlphaPrior[1], "GroupNoiseAlphaPrior[1]", GroupNoiseAlphaPrior[1]);

        parameters.readInto(incHighFreqStoc, "incHighFreqStoc",incHighFreqStoc);
        parameters.readInto(HighFreqStocPrior[0], "HighFreqStocPrior[0]", HighFreqStocPrior[0]);
        parameters.readInto(HighFreqStocPrior[1], "HighFreqStocPrior[1]", HighFreqStocPrior[1]);

        parameters.readInto(incProfileEvo, "incProfileEvo",incProfileEvo);
        parameters.readInto(FitEvoExponent, "FitEvoExponent",FitEvoExponent);
        parameters.readInto(ProfileEvoPrior[0], "ProfileEvoPrior[0]", ProfileEvoPrior[0]);
        parameters.readInto(ProfileEvoPrior[1], "ProfileEvoPrior[1]", ProfileEvoPrior[1]);

        parameters.readInto(incProfileEnergyEvo, "incProfileEnergyEvo",incProfileEnergyEvo);
        parameters.readInto(ProfileEnergyEvoPrior[0], "ProfileEnergyEvoPrior[0]", ProfileEnergyEvoPrior[0]);
        parameters.readInto(ProfileEnergyEvoPrior[1], "ProfileEnergyEvoPrior[1]", ProfileEnergyEvoPrior[1]);

	parameters.readInto(incWideBandNoise, "incWideBandNoise",incWideBandNoise);

        parameters.readInto(incProfileFit, "incProfileFit",incProfileFit);
        parameters.readInto(EvoRefFreq, "EvoRefFreq",EvoRefFreq);
        parameters.readInto(ProfileFitPrior[0], "ProfileFitPrior[0]", ProfileFitPrior[0]);
        parameters.readInto(ProfileFitPrior[1], "ProfileFitPrior[1]", ProfileFitPrior[1]);

        parameters.readInto(FitLinearProfileWidth, "FitLinearProfileWidth",FitLinearProfileWidth);
        parameters.readInto(LinearProfileWidthPrior[0], "LinearProfileWidthPrior[0]", LinearProfileWidthPrior[0]);
        parameters.readInto(LinearProfileWidthPrior[1], "LinearProfileWidthPrior[1]", LinearProfileWidthPrior[1]);

        parameters.readInto(incDMEQUAD, "incDMEQUAD",incDMEQUAD);
        parameters.readInto(DMEQUADPrior[0], "DMEQUADPrior[0]", DMEQUADPrior[0]);
        parameters.readInto(DMEQUADPrior[1], "DMEQUADPrior[1]", DMEQUADPrior[1]);

	
	parameters.readInto(offPulseLevel, "offPulseLevel",offPulseLevel);

        parameters.readInto(strBuf, "ProfFile", string("Prof.dat"));
        strcpy(ProfFile, strBuf.data());

	parameters.readInto(numProfComponents, "numProfComponents",numProfComponents);

        parameters.readInto(incWidthJitter, "incWidthJitter", incWidthJitter);
        parameters.readInto(WidthJitterPrior[0], "WidthJitterPrior[0]", WidthJitterPrior[0]);
        parameters.readInto(WidthJitterPrior[1], "WidthJitterPrior[1]", WidthJitterPrior[1]);

        parameters.readInto(JitterProfComp, "JitterProfComp", JitterProfComp);

        parameters.readInto(debug, "debug", debug);

        parameters.readInto(ProfileBaselineTerms, "ProfileBaselineTerms", ProfileBaselineTerms);


        parameters.readInto(incProfileNoise, "incProfileNoise", incProfileNoise);
        parameters.readInto(ProfileNoiseCoeff, "ProfileNoiseCoeff", ProfileNoiseCoeff);
        parameters.readInto(ProfileNoiseAmpPrior[0], "ProfileNoiseAmpPrior[0]", ProfileNoiseAmpPrior[0]);
        parameters.readInto(ProfileNoiseAmpPrior[1], "ProfileNoiseAmpPrior[1]", ProfileNoiseAmpPrior[1]);
        parameters.readInto(ProfileNoiseSpecPrior[0], "ProfileNoiseSpecPrior[0]", ProfileNoiseSpecPrior[0]);
        parameters.readInto(ProfileNoiseSpecPrior[1], "ProfileNoiseSpecPrior[1]", ProfileNoiseSpecPrior[1]);
	
	parameters.readInto(SubIntToFit, "SubIntToFit", SubIntToFit);
	parameters.readInto(ChannelToFit, "ChannelToFit", ChannelToFit);


	parameters.readInto(NProfileEvoPoly, "NProfileEvoPoly", NProfileEvoPoly);

	parameters.readInto(usecosiprior, "usecosiprior", usecosiprior);
	parameters.readInto(incWidthEvoTime, "incWidthEvoTime", incWidthEvoTime);
	parameters.readInto(incExtraProfComp, "incExtraProfComp", incExtraProfComp);

        parameters.readInto(removeBaseline, "removeBaseline", removeBaseline);


	parameters.readInto(incPrecession, "incPrecession", incPrecession);

	parameters.readInto(incTimeCorrProfileNoise, "incTimeCorrProfileNoise", incTimeCorrProfileNoise);

	parameters.readInto(phasePriorExpansion, "phasePriorExpansion", phasePriorExpansion);

	parameters.readInto(ProfileNoiseMethod, "ProfileNoiseMethod", ProfileNoiseMethod);
	
	parameters.readInto(FitPrecAmps, "FitPrecAmps", FitPrecAmps);

	parameters.readInto(NProfileTimePoly, "NProfileTimePoly", NProfileTimePoly);



	parameters.readInto(incProfileScatter, "incProfileScatter", incProfileScatter);	
	parameters.readInto(ScatterPBF, "ScatterPBF", ScatterPBF);
	
    } catch(ConfigFile::file_not_found oError) {
        printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
    } // try

}

void setTNPriors(char *ConfigFileName, double **Dpriors, long double **TempoPriors, int TPsize, int DPsize){

//This function overwrites the default values for the priors sent to multinest, and the long double priors used by tempo2, you need to be aware of what dimension is what if you use this function.

//THe order of the parameters is always the same:
//Timing Model parameters (linear or non linear)
//Jumps
//EFAC(s) 
//EQUAD
//Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)

	for(int i =0;i<TPsize; i++){	
	//	printf("TP %i \n", i);

	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "TempoPriors[%i][0]", i);
		parameters.readInto(TempoPriors[i][0], buffer, TempoPriors[i][0]);
		n=sprintf (buffer, "TempoPriors[%i][1]", i);
                parameters.readInto(TempoPriors[i][1], buffer, TempoPriors[i][1]);
		n=sprintf (buffer, "TempoPriors[%i][2]", i);
		parameters.readInto(TempoPriors[i][2], buffer, TempoPriors[i][2]);
		n=sprintf (buffer, "TempoPriors[%i][3]", i);
                parameters.readInto(TempoPriors[i][3], buffer, TempoPriors[i][3]);


	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}



	for(int i =0;i<DPsize; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
                char buffer [50];
                int n;
                n=sprintf (buffer, "Dpriors[%i][0]", i);
				parameters.readInto(Dpriors[i][0], buffer, Dpriors[i][0]);
				n=sprintf (buffer, "Dpriors[%i][1]", i);
				parameters.readInto(Dpriors[i][1], buffer, Dpriors[i][1]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}





}


void setFrequencies(char *ConfigFileName, double *SampleFreq, int numRedfreqs, int numDMfreqs, int numRedLogFreqs, int numDMLogFreqs, double RedLowFreq, double DMLowFreq, double RedMidFreq, double DMMidFreq){

//This function sets or overwrites the default values for the sampled frequencies sent to multinest


	int startpoint=0;
	double RedLogDiff = log10(RedMidFreq) - log10(RedLowFreq);
	for(int i =0; i < numRedLogFreqs; i++){
		SampleFreq[startpoint]=pow(10.0, log10(RedLowFreq) + i*RedLogDiff/numRedLogFreqs);
		startpoint++;
		printf("%i %g %g \n", i, log10(RedLowFreq) - i*log10(RedLowFreq)/numRedLogFreqs, SampleFreq[startpoint-1]);
		
	}
	
	for(int i =0;i < numRedfreqs-numRedLogFreqs; i++){	
		SampleFreq[startpoint]=i+RedMidFreq;
		startpoint++;
		//printf("making freqs %i %g\n", startpoint-1, SampleFreq[startpoint-1]);
	
	}
        for(int i =0;i < numDMfreqs; i++){
                SampleFreq[startpoint]=i+1;
		startpoint++;
		//printf("making freqs %i %g", startpoint+i, SampleFreq[startpoint+i]);
        }

	startpoint=0;
	for(int i =0;i<numRedfreqs; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "SampleFreq[%i]", i);
		parameters.readInto(SampleFreq[i], buffer, SampleFreq[i]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	startpoint=startpoint+numRedfreqs;
	for(int i =0;i<numDMfreqs; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "SampleFreq[%i]", startpoint+i);
		parameters.readInto(SampleFreq[startpoint+i], buffer, SampleFreq[startpoint+i]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}
}



void GetGroupsToFit(char *ConfigFileName, int incGroupNoise, int **FitForGroup, int incBandNoise, int **FitForBand){

//This function reads in the groups that will be fit as Group Noise terms



	for(int i =0;i<incGroupNoise; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			/* We can check whether a value is not set in the file by doing
			* if(! parameters.readInto(variable, "name", default)) {
			*   printf("WARNING");
			* }
			*
			* At the moment I was too lazy to print warning messages, and the
			* default value from this file is used in that case.
			*
			* Note: the timing model parameters are not done implemented yet
			*/
			char buffer [50];
			int n;
			n=sprintf (buffer, "FitForGroup[%i][0]", i);
			parameters.readInto(FitForGroup[i][0], buffer, FitForGroup[i][0]);
                        n=sprintf (buffer, "FitForGroup[%i][1]", i);
                        parameters.readInto(FitForGroup[i][1], buffer, FitForGroup[i][1]);
                        n=sprintf (buffer, "FitForGroup[%i][2]", i);
                        parameters.readInto(FitForGroup[i][2], buffer, FitForGroup[i][2]);
                        n=sprintf (buffer, "FitForGroup[%i][3]", i);
                        parameters.readInto(FitForGroup[i][3], buffer, FitForGroup[i][3]);
                        n=sprintf (buffer, "FitForGroup[%i][4]", i);
                        parameters.readInto(FitForGroup[i][4], buffer, FitForGroup[i][4]);
                        n=sprintf (buffer, "FitForGroup[%i][5]", i);
                        parameters.readInto(FitForGroup[i][5], buffer, FitForGroup[i][5]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	for(int i =0;i<incBandNoise; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			/* We can check whether a value is not set in the file by doing
			* if(! parameters.readInto(variable, "name", default)) {
			*   printf("WARNING");
			* }
			*
			* At the moment I was too lazy to print warning messages, and the
			* default value from this file is used in that case.
			*
			* Note: the timing model parameters are not done implemented yet
			*/
			char buffer [50];
			int n;
			n=sprintf (buffer, "FitForBand[%i][0]", i);
			parameters.readInto(FitForBand[i][0], buffer, FitForBand[i][0]);
                        n=sprintf (buffer, "FitForBand[%i][1]", i);
                        parameters.readInto(FitForBand[i][1], buffer, FitForBand[i][1]);
                        n=sprintf (buffer, "FitForBand[%i][2]", i);
                        parameters.readInto(FitForBand[i][2], buffer, FitForBand[i][2]);
                        n=sprintf (buffer, "FitForBand[%i][3]", i);
                        parameters.readInto(FitForBand[i][3], buffer, FitForBand[i][3]);

		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}
}



void GetProfileFitInfo(char *ConfigFileName, int numProfComponents, int *numGPTAshapecoeff, int *numProfileFitCoeff, int *numEvoCoeff, int *numFitEvoCoeff, int *numGPTAstocshapecoeff, double *ProfCompSeps, double &TemplateChanWidth, int *numTimeCorrCoeff, int incExtraComp, double **FitForExtraComp, int *FitCompWidths, int *FitCompPos){

//This function reads in the groups that will be fit as Group Noise terms



	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numGPTAshapecoeff[%i]", i);
			parameters.readInto(numGPTAshapecoeff[i], buffer, numGPTAshapecoeff[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "ProfCompSeps[%i]", i);
			parameters.readInto(ProfCompSeps[i], buffer, ProfCompSeps[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numProfileFitCoeff[%i]", i);
			parameters.readInto(numProfileFitCoeff[i], buffer, numProfileFitCoeff[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numEvoCoeff[%i]", i);
			parameters.readInto(numEvoCoeff[i], buffer, numEvoCoeff[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numFitEvoCoeff[%i]", i);
			parameters.readInto(numFitEvoCoeff[i], buffer, numFitEvoCoeff[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numGPTAstocshapecoeff[%i]", i);
			parameters.readInto(numGPTAstocshapecoeff[i], buffer, numGPTAstocshapecoeff[i]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

	try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		parameters.readInto(TemplateChanWidth, "TemplateChanWidth", TemplateChanWidth);
	} 
	catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	} // try


	for(int i =0;i<numProfComponents; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			char buffer [50];
			int n;
			n=sprintf (buffer, "numTimeCorrCoeff[%i]", i);
			parameters.readInto(numTimeCorrCoeff[i], buffer, numTimeCorrCoeff[i]);

			n=sprintf (buffer, "FitCompWidths[%i]", i);
			parameters.readInto(FitCompWidths[i], buffer, FitCompWidths[i]);

			n=sprintf (buffer, "FitCompPos[%i]", i);
			parameters.readInto(FitCompPos[i], buffer, FitCompPos[i]);



		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	for(int i =0;i<incExtraComp; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = ConfigFileName;//string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			/* We can check whether a value is not set in the file by doing
			* if(! parameters.readInto(variable, "name", default)) {
			*   printf("WARNING");
			* }
			*
			* At the moment I was too lazy to print warning messages, and the
			* default value from this file is used in that case.
			*
			* Note: the timing model parameters are not done implemented yet
			*/
			char buffer [50];
			int n;
			n=sprintf (buffer, "FitForExtraComp[%i][0]", i);
			parameters.readInto(FitForExtraComp[i][0], buffer, FitForExtraComp[i][0]);   //Comp Type
                        n=sprintf (buffer, "FitForExtraComp[%i][1]", i);
                        parameters.readInto(FitForExtraComp[i][1], buffer, FitForExtraComp[i][1]);  //NumCoeff
                        n=sprintf (buffer, "FitForExtraComp[%i][2]", i);
                        parameters.readInto(FitForExtraComp[i][2], buffer, FitForExtraComp[i][2]);  //Include Phase Equivalent
//                        n=sprintf (buffer, "FitForBand[%i][3]", i);
  //                      parameters.readInto(FitForBand[i][3], buffer, FitForBand[i][3]);

		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}

}

void setShapePriors(char *ConfigFileName, double **ShapePriors, double **BetaPrior, int numcoeff, int numcomps){


	for(int i =0;i<numcoeff; i++){	
	//	printf("TP %i \n", i);

	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);
		char buffer [50];
		int n;
		
  		n=sprintf (buffer, "ShapePriors[%i][0]", i);
		parameters.readInto(ShapePriors[i][0], buffer, ShapePriors[i][0]);
		n=sprintf (buffer, "ShapePriors[%i][1]", i);
                parameters.readInto(ShapePriors[i][1], buffer, ShapePriors[i][1]);


	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}



    // Use a configfile, if we can, to overwrite the defaults set in this file.
    for(int i = 0; i < numcomps; i++){
	    try {
		string strBuf;
		strBuf = ConfigFileName;//string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		char buffer [50];
                int n;


		n=sprintf (buffer, "BetaPrior[%i][0]", i);
		parameters.readInto(BetaPrior[i][0], buffer, BetaPrior[i][0]);
		n=sprintf (buffer, "BetaPrior[%i][1]", i);
		parameters.readInto(BetaPrior[i][1], buffer, BetaPrior[i][1]);


	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try
	}





}
