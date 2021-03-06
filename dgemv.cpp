#include <config.h>
#include <math.h>
#include <algorithm>
#include "dgemv.h"
#include "TKlog.h"

#ifdef HAVE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
void vector_dgemv(double *A, double *vecin,double *vecout,int rowa, int cola, char AT)
{

        int M,N,K;

        double alpha=1;
        double beta=0;
        int incX=1;
        int incY=1;

        CBLAS_TRANSPOSE at;

        if(AT == 'N'){
          at = CblasNoTrans;
        }
        else if(AT == 'T'){
          at = CblasTrans;
        }

	logtchk("Running MKL dgemv");
        cblas_dgemv(CblasColMajor, at, rowa, cola, alpha, A, rowa, vecin, incX, beta, vecout, incY);
	logtchk("Finished MKL dgemv");
}
#else
void vector_dgemv(double *A, double *vecin,double *vecout,int rowa, int cola, char AT)
{

	int M,N,K;

	double alpha=1;
	double beta=0;
	int incX=1;
	int incY=1;

	logtchk("Running BLAS dgemv");
	dgemv_(&AT, &rowa, &cola, &alpha, A, &rowa, vecin, &incX, &beta, vecout, &incY);
	logtchk("Finished BLAS dgemv");
}
#endif

void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT)
{

	int M,N,K;
	double *a;

	double alpha=1;
	double beta=0;
	int incX=1;
	int incY=1;

	a = dgemv_ctof(A, rowa, cola); 
	
	dgemv_(&AT, &rowa, &cola, &alpha, a, &rowa, vecin, &incX, &beta, vecout, &incY);
  
  delete a;
}


double* dgemv_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dgemv_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
