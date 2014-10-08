/*
Abtin Rasoulian
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"

void cpd_pr(
		double* x,
		double* y, 
        double* sigma2,
        double* P,
        int N,
		int M,
        int D
        )

{
  int		n, m, d;
  double	ksig, diff, razn, sp;
  
  ksig = -2.0 * *sigma2;
  
  for (n=0; n < N; n++) {
      sp=0;
      for (m=0; m < M; m++) {
          razn=0;
          for (d=0; d < D; d++) {
             diff=*(x+n+d*N)-*(y+m+d*M);  
             diff=diff*diff;
             razn+=diff;
          }
          
          *(P+m+n*M)=exp(razn/ksig);
          sp+=*(P+m+n*M);
      }
      for (m=0; m < M; m++) {
          *(P+m+n*M)=*(P+m+n*M)/ sp;
      }
  }
  return;
}

/* Input arguments */
#define IN_x		prhs[0]
#define IN_y		prhs[1]
#define IN_sigma2	prhs[2]


/* Output arguments */
#define OUT_P		plhs[0]


/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *sigma2, *P;
  int     N, M, D;
  
  /* Get the sizes of each input argument */
  N = mxGetM(IN_x);
  M = mxGetM(IN_y);
  D = mxGetN(IN_x);
  
  /* Create the new arrays and set the output pointers to them */
  OUT_P     = mxCreateDoubleMatrix(M, N, mxREAL);

  /* Assign pointers to the input arguments */
  x      = mxGetPr(IN_x);
  y       = mxGetPr(IN_y);
  sigma2       = mxGetPr(IN_sigma2);

 
  
  /* Assign pointers to the output arguments */
  P      = mxGetPr(OUT_P);

   
  /* Do the actual computations in a subroutine */
  cpd_pr(x, y, sigma2, P, N, M, D);
  
  return;
}


