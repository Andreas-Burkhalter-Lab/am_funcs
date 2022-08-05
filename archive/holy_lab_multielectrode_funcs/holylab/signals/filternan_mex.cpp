#include "mex.h"
#include <math.h>


/*
 * filternan: filter with an IIR along one dimension, ignoring NaNs
 *
 * Syntax:
 *   imout = filternan(b,a,im);
 *   imout = filternan(b,a,im,zi);
 * where
 *   im is the input image (may be multidimensional)
 * and
 *   imout is the filtered image.
 *
 * Copyright 2006 by Timothy E. Holy
 */

void filternan_work(const float *im,int n_dims,const int *sz_im,const double *b, const double *a,const double *zi,double *ziwork,int lfilt,float *imout);
int validate_dimensions(const int *,const int *,int);
int isVector(const int *sz,const int n_dims);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *im;
  const double *b,*a;
  double *zi,*ziwork;
  float *imout;
  int n_dims,dimIndex,lb,la,lzi;
  const int *sz_im;
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("filternan_mex: requires three or four inputs");
  if (nlhs != 1)
    mexErrMsgTxt("filternan_mex: requires one output");

  // Parse the inputs
  // filter b
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("filternan_mex: filter b must be a real double-precision vector");
  b = mxGetPr(curarg);
  if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
    mexErrMsgTxt("filternan_mex: filter b must be a vector");
  lb = mxGetNumberOfElements(curarg);

  // filter a
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("filternan_mex: filter a must be a real double-precision vector");
  a = mxGetPr(curarg);
  if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
    mexErrMsgTxt("filternan_mex: filter a must be a vector");
  la = mxGetNumberOfElements(curarg);

  // image
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("filternan_mex: data must be a real single-precision array");
  im = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz_im = mxGetDimensions(curarg);

  // filter zi
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
      mexErrMsgTxt("filternan_mex: filter zi must be a real double-precision vector");
    zi = mxGetPr(curarg);
    if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
      mexErrMsgTxt("filternan_mex: filter zi must be a vector");
    lzi = *mxGetDimensions(curarg);
  }
  else {
    // Initialize with zeros
    lzi = la - 1;
    zi = (double*) mxCalloc(lzi,sizeof(double));
  }

  ziwork = (double*) mxMalloc(lzi*sizeof(double));

  // Size checking
  if (la != lb || la != lzi+1)
    mexErrMsgTxt("filternan_mex: lengths of b, a, and zi must match");

  // Set up the output
  plhs[0] = mxCreateNumericArray(n_dims,sz_im,mxSINGLE_CLASS,mxREAL);
  imout = (float *) mxGetData(plhs[0]);

  // Do the actual work
  filternan_work(im,n_dims,sz_im,b,a,zi,ziwork,lb,imout);

  if (nrhs == 3)
    mxFree(zi);
  mxFree(ziwork);

  return;
}


/*
 * Do the filtering.
 */
void filternan_work(const float *im,int n_dims,const int *sz_im,const double *b, const double *a,const double *zi,double *ziwork,int lfilt,float *imout)
{
  int n_pixels,colIterator,ziIterator,dimIterator,initialize;
  double *ziworktmp;
  const float *im_end;
  const double *zi_end, *zitmp;
  double tmp;

  // Compute the total number of pixels in the image
  n_pixels = 1;
  dimIterator = 0;
  while (dimIterator < n_dims)
    n_pixels *= sz_im[dimIterator++];
  im_end = im+n_pixels;

  zi_end = zi+lfilt-1;

  //  mexPrintf("zi_end-zi %d, lfilt %d\n",zi_end-zi,lfilt);

  for (; im < im_end; ) {
    initialize = 1;
    for (colIterator = 0; colIterator < *sz_im; colIterator++, im++, imout++) {
      if (!mxIsNaN(*im)) {
	if (initialize) {
	  initialize = 0;
	  for (zitmp = zi, ziworktmp = ziwork; zitmp < zi_end; zitmp++, ziworktmp++) {
	    *ziworktmp = *zitmp * *im;
	  }
	}
	tmp = b[0] * *im + ziwork[0];
	for (ziIterator = 1; ziIterator < lfilt-1; ziIterator++)
	  ziwork[ziIterator-1] = b[ziIterator] * *im + ziwork[ziIterator] -
	    a[ziIterator] * tmp;
	ziwork[lfilt-2] = b[lfilt-1] * *im - a[lfilt-1] * tmp;
	*imout = tmp;
      } else {
	initialize = 1;
	*imout = (float) mxGetNaN();
      }
    }  // end colIterator
  } // end loop over columns
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
