/* IIR gaussian filtering. This is slightly modified from code
   provided by Christoph Lampert (see copyright info below), based on
   the Triggs_Sdika paper. */
/* Main modifications:
 *   1. Matlab wrapper
 *   2. Have it work on columns of images (compute coefficients only once)
 *   3. NaN-ifying (ignore NaNs at edges---this is done in the parent function)
 *   4. Set boundary conditions to 0; then one can get more edge
        reasonable behavior by also doing filtering on an image of 1s
        to compute a "normalizing" factor (see imfilter_gaussian).
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
//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "imiterators.cxx"

#define SIGMA_IIR_THRESH 2

//typedef float real;

// forward declarations
void imfg_work_single(float *imout,const int n_dims,const int *sz_im,const double *sigma);
void imfg_work_double(double *imout,const int n_dims,const int *sz_im,const double *sigma);
template <class real>
static void FilterCoefficients_and_TriggsCoefficients(real *b, double *M, double sigma);
template <class dataType>
void imfg_fir(colIterator<dataType> ci,dataType *tmp,dataType *h,int hlen);
template <class dataType>
void imfg_iir(colIterator<dataType> ci,const dataType b[],const double M[]);
template <class dataType>
void imfg_work(dataType *imout,const int n_dims,const int *sz_im,const double *sigma);
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
  const double *sigma;
  void *imout;
  int n_dims;
  size_t n_pixels;
  const int *sz_im;
  const mxArray *curarg;

  if (nrhs != 2)
    mexErrMsgTxt("imfilter_gaussian_mex: requires two inputs");
  if (nlhs != 1)
    mexErrMsgTxt("imfilter_gaussian_mex: requires one output");

  // Parse the inputs
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("imfilter_gaussian_mex: data must be a real array");
  im = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz_im = mxGetDimensions(curarg);
  n_pixels = mxGetNumberOfElements(curarg);

  // sigma
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("imfilter_gaussian_mex: sigma must be a real double-precision vector");
  if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
    mexErrMsgTxt("imfilter_gaussian_mex: sigma must be a vector");
  if (mxGetNumberOfElements(curarg) != n_dims)
    mexErrMsgTxt("imfilter_gaussian_mex: number of elements in sigma must match image dimensionality");
  sigma = mxGetPr(curarg);

  // Set up the output
  plhs[0] = mxCreateNumericArray(n_dims,sz_im,mxGetClassID(prhs[0]),mxREAL);
  imout = mxGetData(plhs[0]);

  // Copy over the image data (filtering is done in-place)
  memcpy(imout,im,n_pixels*mxGetElementSize(prhs[0]));

  // Do the actual work
  if (mxIsDouble(prhs[0]))
    imfg_work_double((double *)imout,n_dims,sz_im,sigma);
  else
    imfg_work_single((float *)imout,n_dims,sz_im,sigma);

  return;
}

template <class dataType>
void imfg_work(dataType *imout,const int n_dims,const int *sz_im,const double *sigma)
{
  int dimIndex,hlen,i;
  bool use_iir;
  dataType *tmp,*h;
  dataType hsum;

  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    if (sigma[dimIndex] > 0 && sz_im[dimIndex] > 1) {
      multiDim1dIterator<dataType> mdi(imout,sz_im,n_dims,dimIndex);  // optional range_start & range_end for multithreading
      use_iir = sigma[dimIndex] > SIGMA_IIR_THRESH && sz_im[dimIndex] > 3;
      if (use_iir) {
	// IIR filtering. This uses the Triggs method.
	dataType b[6];
	double M[9];
	FilterCoefficients_and_TriggsCoefficients(b,M,sigma[dimIndex]);
	for (; !mdi.at_end(); mdi.inc_col()) {
	  imfg_iir(mdi.col_iterator(),b,M);
	}
      } else {
	// FIR filtering. This is just the direct method
	// Create temporary space, so that filtered values do not get used prematurely
	tmp = (dataType *) mxMalloc(sz_im[dimIndex]*sizeof(dataType));
	// Create the gaussian filter
	hlen = (int) ceil(4*sigma[dimIndex]);
	if (hlen > sz_im[dimIndex])
	  hlen = sz_im[dimIndex];
	h = (dataType *) mxMalloc((2*hlen+1)*sizeof(dataType));
	h += hlen;
	hsum = 0;
	for (i = 0; i <= hlen; i++) {
	  h[i] = exp(-((dataType)(i*i))/(2.0*sigma[dimIndex]*sigma[dimIndex]));
	  h[-i] = h[i];
	  if (i > 0)
	    hsum += 2*h[i];
	  else
	    hsum += h[i];
	}
	for (i = 0; i <= hlen; i++) {
	  h[i] /= hsum;
	  if (i > 0)
	    h[-i] /= hsum;
	}
	for (; !mdi.at_end(); mdi.inc_col()) {
	  imfg_fir(mdi.col_iterator(),tmp,h,hlen);
	}
	mxFree(h-hlen);
	mxFree(tmp);
      }
    }
  }
}

void imfg_work_single(float *imout,const int n_dims,const int *sz_im,const double *sigma)
{
  imfg_work<float>(imout,n_dims,sz_im,sigma);
}

void imfg_work_double(double *imout,const int n_dims,const int *sz_im,const double *sigma)
{
  imfg_work<double>(imout,n_dims,sz_im,sigma);
}



template <class real>
void imfg_iir(colIterator<real> ci,const real b[],const double M[])
{
  real iminus;    //  hypothetical signal value before start time, i.e. t=-1,-2,-3 
  real iplus;	  // last sample of original signal, needed for right hand boundary values

  real v[3],uu[3];
  double p0,p1,p2,p3;

  /*
  colIterator<float> citmp = ci;
  for (citmp; citmp < citmp.end(); citmp++)
    mexPrintf("%g ",*citmp);
  mexPrintf("\n");
  colIterator<float> citmp2 = ci;
  */

  //iplus = ci[ci.length()-1];  // Last sample of original signal
  iplus = 0;
  
  // Initialize first three points in the causal direction
  //iminus = *ci;
  iminus = 0;
  //mexPrintf("iplus %g, iminus %g\n",iplus,iminus);
  *ci -= -b[1]*(iminus-*ci) + b[2]*(iminus-*ci) - b[3]*(iminus-*ci);
  ci++;
  *ci -= -b[1]*(ci[-1]-*ci) + b[2]*(iminus-*ci) - b[3]*(iminus-*ci); // second element
  ci++;
  *ci -= -b[1]*(ci[-1]-*ci) + b[2]*(ci[-2]-*ci) - b[3]*(iminus-*ci); // third element
  ci++;

  // the remainder of the causal filter
  for (p0 = *ci,p1 = ci[-1],p2 = ci[-2],p3 = ci[-3]; ci < ci.end(); p3=p2,p2=p1,p1=p0,ci++) {
    // Do this in a way that preserves double-precision accuracy
    p0 = *ci;
    p0 -= -b[1]*(p1-p0) + b[2]*(p2-p0) - b[3]*(p3-p0);
    *ci = (float) p0;
  }
  // note ci is now positioned at the beyond-the-end position, so
  // ci[-1] corresponds to the last element
  /*
  for (; citmp2 < citmp2.end(); citmp2++)
    mexPrintf("%g ",*citmp2);
  mexPrintf("\n");
  */

  // The first three steps of the anticausal filter 
  uu[0] = ci[-1]-iplus;
  uu[1] = ci[-2]-iplus;
  uu[2] = ci[-3]-iplus;
  v[0] = M[0]*uu[0] - M[1]*uu[1] + M[2]*uu[2] + iplus;
  v[1] = M[3]*uu[0] - M[4]*uu[1] + M[5]*uu[2] + iplus;
  v[2] = M[6]*uu[0] - M[7]*uu[1] + M[8]*uu[2] + iplus;
  //mexPrintf("v: %g %g %g\n",v[0],v[1],v[2]);

  ci--;
  *ci -= -b[1]*(v[0]-*ci) + b[2]*(v[1]-*ci) - b[3]*(v[2]-*ci);
  p3 = *ci;
  ci--;
  *ci -= -b[1]*(ci[1]-*ci) + b[2]*(v[0]-*ci) - b[3]*(v[1]-*ci);
  p2 = *ci;
  ci--;
  *ci -= -b[1]*(ci[1]-*ci) + b[2]*(ci[2]-*ci) - b[3]*(v[0]-*ci);
  p1 = *ci;
  ci--;
  //mexPrintf("%g %g %g\n",p1,p2,p3);

  // The rest of the anticausal filter
  for (p0=*ci,p1=ci[1],p2=ci[2],p3=ci[3]; ci >= ci.begin(); p3=p2,p2=p1,p1=p0,ci--) {
    p0 = *ci;
    p0 -= -b[1]*(p1-p0) + b[2]*(p2-p0) - b[3]*(p3-p0);
    *ci = (real) p0;
  }
}



template <class dataType>
void imfg_fir(colIterator<dataType> ci,dataType *tmp,dataType *h,int hlen)
{
  int colLength,pixCounter,left,right;
  colIterator<dataType> pixIterator;
  dataType *tmptmp,*htmp,*htmpend;
  
  colLength = ci.length();
  for (pixCounter = 0,tmptmp = tmp; ci < ci.end(); ci++,pixCounter++,tmptmp++) {
    *tmptmp = 0;
    // Compute the margins on either side, to make sure we don't use invalid data
    left = (pixCounter < hlen) ? pixCounter : hlen;
    right = colLength - pixCounter - 1;
    right = (right < hlen) ? right : hlen;
    // Compute the gaussian filtered data
    for (pixIterator = ci - left, htmp = h-left, htmpend = h+right; htmp <= htmpend; pixIterator++,htmp++)
      *tmptmp += *pixIterator * *htmp;
  }
  // Now copy the result back to the image data
  for (ci = ci.begin(); ci < ci.end(); ci++,tmp++) {
    *ci = *tmp;
  }
}


template <class real>
static void FilterCoefficients_and_TriggsCoefficients(real *b, double *M, double sigma)
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


int isVector(const int *sz,const int n_dims)
{
  int n_nonunity,dimIndex;

  n_nonunity = 0;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (sz[dimIndex] > 1)
      n_nonunity++;
  return (n_nonunity < 2);
}
