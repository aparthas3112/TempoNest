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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "configfile.h"

void setupparams(char* ConfigFileName, char* root, int& numTempo2its, int& incEFAC, int& incEQUAD,
                 int& incRED, int& incDM, double* EFACPrior, double* EQUADPrior, double* AlphaPrior,
                 double* AmpPrior, double* DMAlphaPrior, double* DMAmpPrior, double& numRedCoeff,
                 double& numDMCoeff, double& FourierSig, char* whiteflag, int& RedPriorType,
                 int& DMPriorType, int& EQUADPriorType, int& EFACPriorType, int& useOriginalErrors,
                 int& StoreTMatrix, int& debug)
{

    // General parameters:
    // Use GPUs 0=No, 1=Yes

    debug = 0;

    StoreTMatrix = 1;  // Recompute TMatrices when computing new bats - default is just precompute
                       // and use those

    // Root of the results files,relative to the directory in which TempoNest is run. This will be
    // followed by the pulsar name, and then the individual output file extensions.
    strcpy(root, "results/Example1-");
    // numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.
    // Should only be set to 0 if all the priors are set in setTNPriors
    numTempo2its = 1;

    // ModelChoices

    // White noise

    useOriginalErrors = 0;  // Use tempo2 errors before modification by TNEF/TNEQ
    incEFAC =
        0;  // include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
    incEQUAD = 0;  // include EQUAD: 0 = no, 1 = yes

    incRED = 0;  // include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model
                 // independent (L2013)
    incDM = 0;   // include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model
                 // independent (L2013)

    // Priors

    // Remaining priors for the stochastic parameters.

    RedPriorType = 0;    // 0 = Log, 1 = Uniform
    DMPriorType = 0;     // 0 = Log, 1 = Uniform
    EQUADPriorType = 0;  // 0 = Log, 1 = Uniform
    EFACPriorType = 0;   // 0 = Log, 1 = Uniform

    EFACPrior[0] = 0.1;
    EFACPrior[1] = 10;

    EQUADPrior[0] = -10;
    EQUADPrior[1] = -5;

    numRedCoeff = 10;
    numDMCoeff = 10;

    AlphaPrior[0] = 1.1;
    AlphaPrior[1] = 6.1;

    AmpPrior[0] = -20;
    AmpPrior[1] = -10;

    DMAlphaPrior[0] = 1.1;
    DMAlphaPrior[1] = 6.1;

    DMAmpPrior[0] = -18;
    DMAmpPrior[1] = -8;

    FourierSig = 5;

    // Use a configfile, if we can, to overwrite the defaults set in this file.
    try {
        string strBuf;

        strBuf = ConfigFileName;  // string("defaultparameters.conf");
        // string NewstrBuf = ConfigFileName;

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
        parameters.readInto(StoreTMatrix, "StoreTMatrix", StoreTMatrix);

        parameters.readInto(strBuf, "root", string("results/Example1"));
        strcpy(root, strBuf.data());
        parameters.readInto(numTempo2its, "numTempo2its", numTempo2its);
        parameters.readInto(incEFAC, "incEFAC", incEFAC);
        parameters.readInto(incEQUAD, "incEQUAD", incEQUAD);

        parameters.readInto(incRED, "incRED", incRED);
        parameters.readInto(incDM, "incDM", incDM);
        parameters.readInto(EFACPrior[0], "EFACPrior[0]", EFACPrior[0]);
        parameters.readInto(EFACPrior[1], "EFACPrior[1]", EFACPrior[1]);
        parameters.readInto(EQUADPrior[0], "EQUADPrior[0]", EQUADPrior[0]);
        parameters.readInto(EQUADPrior[1], "EQUADPrior[1]", EQUADPrior[1]);
        parameters.readInto(AlphaPrior[0], "AlphaPrior[0]", AlphaPrior[0]);
        parameters.readInto(AlphaPrior[1], "AlphaPrior[1]", AlphaPrior[1]);
        parameters.readInto(AmpPrior[0], "AmpPrior[0]", AmpPrior[0]);
        parameters.readInto(AmpPrior[1], "AmpPrior[1]", AmpPrior[1]);
        parameters.readInto(numRedCoeff, "numRedCoeff", numRedCoeff);
        parameters.readInto(numDMCoeff, "numDMCoeff", numDMCoeff);

        parameters.readInto(DMAlphaPrior[0], "DMAlphaPrior[0]", DMAlphaPrior[0]);
        parameters.readInto(DMAlphaPrior[1], "DMAlphaPrior[1]", DMAlphaPrior[1]);
        parameters.readInto(DMAmpPrior[0], "DMAmpPrior[0]", DMAmpPrior[0]);
        parameters.readInto(DMAmpPrior[1], "DMAmpPrior[1]", DMAmpPrior[1]);

        parameters.readInto(strBuf, "whiteflag", string("-sys"));
        strcpy(whiteflag, strBuf.data());

        parameters.readInto(RedPriorType, "RedPriorType", RedPriorType);
        parameters.readInto(DMPriorType, "DMPriorType", DMPriorType);
        parameters.readInto(EFACPriorType, "EFACPriorType", EFACPriorType);
        parameters.readInto(EQUADPriorType, "EQUADPriorType", EQUADPriorType);
        parameters.readInto(useOriginalErrors, "useOriginalErrors", useOriginalErrors);

        parameters.readInto(debug, "debug", debug);

    } catch (ConfigFile::file_not_found oError) {
        printf("WARNING: parameters file '%s' not found. Using defaults.\n",
               oError.filename.c_str());
    }  // try
}

void setTNPriors(char* ConfigFileName, double** Dpriors, long double** TempoPriors, int TPsize,
                 int DPsize)
{

    // This function overwrites the default values for the priors sent to multinest, and the long
    // double priors used by tempo2, you need to be aware of what dimension is what if you use this
    // function.

    // THe order of the parameters is always the same:
    // Timing Model parameters (linear or non linear)
    // Jumps
    // EFAC(s)
    // EQUAD
    // Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)

    for (int i = 0; i < TPsize; i++) {
        //	printf("TP %i \n", i);

        // Use a configfile, if we can, to overwrite the defaults set in this file.
        try {
            string strBuf;
            strBuf = ConfigFileName;  // string("defaultparameters.conf");
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
            char buffer[50];
            int n;
            n = sprintf(buffer, "TempoPriors[%i][0]", i);
            parameters.readInto(TempoPriors[i][0], buffer, TempoPriors[i][0]);
            n = sprintf(buffer, "TempoPriors[%i][1]", i);
            parameters.readInto(TempoPriors[i][1], buffer, TempoPriors[i][1]);
            n = sprintf(buffer, "TempoPriors[%i][2]", i);
            parameters.readInto(TempoPriors[i][2], buffer, TempoPriors[i][2]);
            n = sprintf(buffer, "TempoPriors[%i][3]", i);
            parameters.readInto(TempoPriors[i][3], buffer, TempoPriors[i][3]);

        } catch (ConfigFile::file_not_found oError) {
            printf("WARNING: parameters file '%s' not found. Using defaults.\n",
                   oError.filename.c_str());
        }  // try
    }

    for (int i = 0; i < DPsize; i++) {

        // Use a configfile, if we can, to overwrite the defaults set in this file.
        try {
            string strBuf;
            strBuf = ConfigFileName;  // string("defaultparameters.conf");
            ConfigFile parameters(strBuf);

            /* We can check whether a value is not set in the file by doing
             * if(! parameters.readInto(variable, "name", default)) {
             *   printf("WARNING");
             * }
             *
             *
             * Note: the timing model parameters are not done implemented yet
             */
            char buffer[50];
            int n;
            n = sprintf(buffer, "Dpriors[%i][0]", i);
            parameters.readInto(Dpriors[i][0], buffer, Dpriors[i][0]);
            n = sprintf(buffer, "Dpriors[%i][1]", i);
            parameters.readInto(Dpriors[i][1], buffer, Dpriors[i][1]);

        } catch (ConfigFile::file_not_found oError) {
            printf("WARNING: parameters file '%s' not found. Using defaults.\n",
                   oError.filename.c_str());
        }  // try
    }
}

void setFrequencies(char* ConfigFileName, double* SampleFreq, int numRedfreqs, int numDMfreqs,
                    int numRedLogFreqs, int numDMLogFreqs, int numScatLogFreqs, double RedLowFreq,
                    double DMLowFreq, double ScatLowFreq, double RedMidFreq, double DMMidFreq,
                    double ScatMidFreq)
{

    // This function sets or overwrites the default values for the sampled frequencies sent to
    // multinest

    int startpoint = 0;
    double RedLogDiff = log10(RedMidFreq) - log10(RedLowFreq);
    for (int i = 0; i < numRedLogFreqs; i++) {
        SampleFreq[startpoint] = pow(10.0, log10(RedLowFreq) + i * RedLogDiff / numRedLogFreqs);
        startpoint++;
        printf("%i %g %g \n", i, log10(RedLowFreq) - i * log10(RedLowFreq) / numRedLogFreqs,
               SampleFreq[startpoint - 1]);
    }

    for (int i = 0; i < numRedfreqs - numRedLogFreqs; i++) {
        SampleFreq[startpoint] = i + RedMidFreq;
        startpoint++;
        // printf("making freqs %i %g\n", startpoint-1, SampleFreq[startpoint-1]);
    }
    for (int i = 0; i < numDMfreqs; i++) {
        SampleFreq[startpoint] = i + 1;
        startpoint++;
        // printf("making freqs %i %g", startpoint+i, SampleFreq[startpoint+i]);
    }

    startpoint = 0;
    for (int i = 0; i < numRedfreqs; i++) {

        // Use a configfile, if we can, to overwrite the defaults set in this file.
        try {
            string strBuf;
            strBuf = ConfigFileName;  // string("defaultparameters.conf");
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
            char buffer[50];
            int n;
            n = sprintf(buffer, "SampleFreq[%i]", i);
            parameters.readInto(SampleFreq[i], buffer, SampleFreq[i]);

        } catch (ConfigFile::file_not_found oError) {
            printf("WARNING: parameters file '%s' not found. Using defaults.\n",
                   oError.filename.c_str());
        }  // try
    }

    startpoint = startpoint + numRedfreqs;
    for (int i = 0; i < numDMfreqs; i++) {

        // Use a configfile, if we can, to overwrite the defaults set in this file.
        try {
            string strBuf;
            strBuf = ConfigFileName;  // string("defaultparameters.conf");
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
            char buffer[50];
            int n;
            n = sprintf(buffer, "SampleFreq[%i]", startpoint + i);
            parameters.readInto(SampleFreq[startpoint + i], buffer, SampleFreq[startpoint + i]);

        } catch (ConfigFile::file_not_found oError) {
            printf("WARNING: parameters file '%s' not found. Using defaults.\n",
                   oError.filename.c_str());
        }  // try
    }
}