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
#include "eigen_config.h"
#include "tempo2.h"

void TNtextOutput(pulsar* psr, int npsr, int newpar, void* context, int ndims,
                  std::vector<double> paramlist, double Evidence, std::string longname,
                  double** paramarray);

// non linear timing model likelihood functions
double NewLRedMarginLogLike(double Cube[], int ndim, double phi[], int nDerived, void* context);

void LRedLikeMNWrap(double* Cube, int& ndim, int& npars, double& lnew, void* context);

void StoreTMatrix();
void getArraySizeInfo();
void OutputMLFiles(int nParameters, double* pdParameterEstimates, double MLike, int startDim);

void readsummary(pulsar* psr, std::string longname, int ndim, void* context, int ndims);
