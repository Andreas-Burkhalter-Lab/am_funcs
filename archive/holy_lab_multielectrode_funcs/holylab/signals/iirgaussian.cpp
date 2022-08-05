/* IIR gaussian filtering. This is only slightly modified from code provided by Christoph Lampert (see copyright info below), based on the Triggs_Sdika paper. */
/* Main modifications:
 *   1. Matlab wrapper
 *   2. Have it work on columns of images (compute coefficients only once)
 *   3. NaN-ifying (ignore NaNs at edges)
 *
 * Modifications copyright 2006 by Timothy E. Holy 
 */



/** recursive_gaussian-float.c
 *  A floating point implementation of 1D recursive Gaussian filtering.
 *
 *  Copyright 2005 Christoph Lampert 
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
**/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

typedef float real;

//#define MAIN
#define MATLABWRAPPER

// forward declarations
static void FilterCoefficients_and_TriggsCoefficients(real *b, double *M, const real sigma);
static void ForwardCoefficients(real *v, const double *M, const real *b, const real *causal, const real iplus, const int len);
void gaussian_filter1d(const real *in, real *out, int len, const real *b, const double *M);
void iirg_work(const float *im,int len,double sigma,float *imout);
int isVector(const int *sz,const int n_dims);

/** Recursive Young-van Vliet Gaussian filter (floating point version)
 * @param in [IN] Signal
 * @param out [OUT] Filtered signal result
 * @parma len [IN] Length of \a in and \a out
 * @param sigma [IN] Standard deviation
 **/
 
/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *im;
  double sigma;
  float *imout;
  int len;
  const mxArray *curarg;

  if (nrhs != 2)
    mexErrMsgTxt("iirgaussian: requires two inputs");
  if (nlhs != 1)
    mexErrMsgTxt("iirgaussian: requires one output");

  // Parse the inputs
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("iirgaussian: data must be a real single-precision array");
  im = (float *) mxGetData(curarg);
  if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
    mexErrMsgTxt("iirgaussian: data must be a vector");
  len = mxGetNumberOfElements(curarg);

  // sigma
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("iirgaussian: sigma must be a real double-precision vector");
  sigma = mxGetScalar(curarg);

  // Set up the output
  plhs[0] = mxCreateNumericMatrix(1,len,mxSINGLE_CLASS,mxREAL);
  imout = (float *) mxGetData(plhs[0]);

  // Do the actual work
  iirg_work(im,len,sigma,imout);

  return;
}

void iirg_work(const float *im,int len,double sigma,float *imout)
{
  real b[6];
  double M[9];

  FilterCoefficients_and_TriggsCoefficients(b, M, sigma);
  gaussian_filter1d(im,imout,len,b,M);
}
  

void gaussian_filter1d(const real *in, real *out, int len, const real *b, const double *M)
{
  int i;
  real iminus;    //  hypothetical signal value before start time, i.e. t=-1,-2,-3 
  real iplus;	  // last sample of original signal, needed for right hand boundary values

  real v[4];

  iplus = in[len-1];
  
  // Initialize first three steps of the causal filter
  iminus = in[0];
  out[0] = in[0] - ( - b[1]*(iminus-in[0]) + b[2]*(iminus-in[0]) - b[3]*(iminus-in[0]) );// (20)
  out[1] = in[1] - ( - b[1]*(out[0]-in[1]) + b[2]*(iminus-in[1]) - b[3]*(iminus-in[1]) );
  out[2] = in[2] - ( - b[1]*(out[1]-in[2]) + b[2]*(out[0]-in[2]) - b[3]*(iminus-in[2]) );
   
  // the remainder of the causal filter
  for (i = 3; i  < len; i++)
	out[i] = in[i] - ( - b[1]*(out[i-1]-in[i]) + b[2]*(out[i-2]-in[i]) - b[3]*(out[i-3]-in[i]) );

  // The first three steps of the anticausal filter 
  ForwardCoefficients(v, M, b, out, iplus, len);
  
  out[len-1] -= - b[1]*(v[0]      -out[len-1]) + b[2]*(v[1]-out[len-1]      ) - b[3]*(v[2]-out[len-1]);
  out[len-2] -= - b[1]*(out[len-1]-out[len-2]) + b[2]*(v[0]-out[len-2]      ) - b[3]*(v[1]-out[len-2]);
  out[len-3] -= - b[1]*(out[len-2]-out[len-3]) + b[2]*(out[len-1]-out[len-3]) - b[3]*(v[0]-out[len-3]);

  // the remainder of the anti-causal filter 
  for (i = len - 4; i >= 0; i--)
	out[i] -= - b[1]*(out[i+1]-out[i]) + b[2]*(out[i+2]-out[i]) - b[3]*(out[i+3]-out[i]);
}

static void FilterCoefficients_and_TriggsCoefficients(real *b, double *M, real sigma)
{
  /** Compute the coefficients for the recursive filter.
   *  @see Young and Van Vliet, Signal Processing (2002)
  **/
  const double m0 = 1.16680; 
  const double m1 = 1.10783;
  const double m2 = 1.40586;
  double b1,b2,b3;
  double q,scale;
  
//  q = 1.31564 * (sqrt(1.+0.490811*sigma*sigma)-1.);	
  if (sigma < 3.3556)
  	q = -0.2568 + 0.5784*sigma + 0.0561*sigma*sigma;
  else
  	q = 2.5091 + 0.9804*(sigma - 3.556);

  scale = (m0+q)*(m1*m1+m2*m2+2*m1*q+q*q);
  b1 = q*(2*m0*m1+m1*m1+m2*m2+(2*m0+4*m1)*q+3*q*q)/scale;
  b2 = q*q*(m0+2*m1+3*q)/scale;
  b3 = q*q*q/scale;

  b[0] = 1.;
  b[1] = b1;
  b[2] = b2;
  b[3] = b3;
  
  /** Compute the matrix M necessary to apply the Triggs boundary conditions.
  *  @see Triggs and Sdika, IEEE Trans. Signal Processing (to appear)
  **/  

  const double a1 = b1; //Triggs paper has different sign conventions 
  const double a2 = -b2; 
  const double a3 = b3; 
  const double norm = 1.0/((1.0+a1-a2+a3)*(1.0+a2+(a1-a3)*a3));
  
  M[0] = norm*(-a3*a1+1.0-a3*a3-a2) ;
  M[1] = -norm*(a3+a1)*(a2+a3*a1) ;
  M[2] = norm*a3*(a1+a3*a2) ;
  M[3] = norm*(a1+a3*a2) ;
  M[4] = norm*(a2-1.0)*(a2+a3*a1 ) ;
  M[5] = -norm*a3*(a3*a1+a3*a3+a2-1.0) ;
  M[6] = norm*(a3*a1+a2+a1*a1-a2*a2) ;
  M[7] = -norm*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3) ;
  M[8] = norm*a3*(a1+a3*a2) ;
}

/** Compute the first three coefficients for the anticausal filter part.
 *  @see Triggs and Sdika
 **/  
static void ForwardCoefficients(real *v, const double *M, const real *b, const real *causal, const real iplus, const int len)
{
  real uplus, vplus;
  real uu[3];
  
  uplus = iplus;
  vplus = uplus;
  
  uu[0] = causal[len-1] - uplus;
  uu[1] = causal[len-2] - uplus;
  uu[2] = causal[len-3] - uplus;
  
  v[0] = M[0]*uu[0] - M[1]*uu[1] + M[2]*uu[2] + vplus;
  v[1] = M[3]*uu[0] - M[4]*uu[1] + M[5]*uu[2] + vplus;
  v[2] = M[6]*uu[0] - M[7]*uu[1] + M[8]*uu[2] + vplus;
}



#ifdef MAIN
int main(int argc, char *argv[])
{
  const int N = 1000;
  int i;
  real signal[N], result[N];
  real sigma;
  real b[6];
  double M[9];

  if (argc != 2)
    return(-1);
  sigma = atof(argv[1]);

  FilterCoefficients_and_TriggsCoefficients(b, M, sigma);

  for (i = 0; i < N; i++)
    signal[i] = 0;   
  //signal[N/2] = 1.;
  signal[10] = 1.;

  gaussian_filter1d(signal, result, N, b, M);
  for (i = 0; i < N; i++)
  {
    printf("%f %f\n", signal[i], result[i]);
  }
 
  return(0);
}
#endif

int isVector(const int *sz,const int n_dims)
{
  int n_nonunity,dimIndex;

  n_nonunity = 0;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (sz[dimIndex] > 1)
      n_nonunity++;
  return (n_nonunity < 2);
}
