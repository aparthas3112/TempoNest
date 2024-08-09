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

void setupparams(char* ConfigFileName, char* root, int& numTempo2its, double& numRedCoeff,
                 double& numDMCoeff, char* whiteflag, int& useOriginalErrors, int& debug)
{

    // General parameters:
    // Use GPUs 0=No, 1=Yes

    debug = 0;

    // Root of the results files,relative to the directory in which TempoNest is run. This will be
    // followed by the pulsar name, and then the individual output file extensions.
    strcpy(root, "results/Example1-");
    // numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.
    // Should only be set to 0 if all the priors are set in setTNPriors
    numTempo2its = 1;

    // ModelChoices

    // White noise

    useOriginalErrors = 0;  // Use tempo2 errors before modification by TNEF/TNEQ

    // Priors

    // Remaining priors for the stochastic parameters.

    numRedCoeff = 10;
    numDMCoeff = 10;

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

        parameters.readInto(strBuf, "root", string("results/Example1"));
        strcpy(root, strBuf.data());
        parameters.readInto(numTempo2its, "numTempo2its", numTempo2its);

        parameters.readInto(numRedCoeff, "numRedCoeff", numRedCoeff);
        parameters.readInto(numDMCoeff, "numDMCoeff", numDMCoeff);

        parameters.readInto(strBuf, "whiteflag", string("-sys"));
        strcpy(whiteflag, strBuf.data());

        parameters.readInto(useOriginalErrors, "useOriginalErrors", useOriginalErrors);

        parameters.readInto(debug, "debug", debug);

    } catch (ConfigFile::file_not_found oError) {
        printf("WARNING: parameters file '%s' not found. Using defaults.\n",
               oError.filename.c_str());
    }  // try
}

void setFrequencies(char* ConfigFileName, double* SampleFreq, int numRedfreqs, int numDMfreqs,
                    int numRedLogFreqs, int numDMLogFreqs, int numScatLogFreqs, double RedLowFreq,
                    double DMLowFreq, double ScatLowFreq, double RedMidFreq, double DMMidFreq,
                    double ScatMidFreq)
{

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
    }
    for (int i = 0; i < numDMfreqs; i++) {
        SampleFreq[startpoint] = i + 1;
        startpoint++;
    }
}