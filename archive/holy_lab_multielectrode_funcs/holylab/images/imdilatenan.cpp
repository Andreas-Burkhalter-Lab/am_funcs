#include "mex.h"
#include <math.h>
#include <string.h>

/*
 * imdilatenan: replace NaN's with neighboring finite values
 *
 * Syntax:
 *   imout = imdilatenan(im)
 * where
 *   im is the input image (may be multidimensional)
 * and
 *   imout is the output image.
 *
 * Copyright 2006 by Timothy E. Holy
 */

void imdilatenan_work(const float *im,int n_dims,const int *sz_im,float *imout);
/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *im;
  const double *h;
  float *imout;//,fillval;
  int n_dims,dimIndex,lh,stepsize;
  const int *sz_im,*sz_h;
  int *coords;
  const mxArray *curarg;

  if (nrhs != 1)
    mexErrMsgTxt("imdilatenan: requires one input");
  if (nlhs != 1)
    mexErrMsgTxt("imdilatenan: requires one output");

  // Parse the inputs
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("imdilatenan: image must be a real single-precision array");
  im = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz_im = mxGetDimensions(curarg);

  // Set up the output
  plhs[0] = mxCreateNumericArray(n_dims,sz_im,mxSINGLE_CLASS,mxREAL);
  imout = (float *) mxGetData(plhs[0]);

  // Do the actual work
  imdilatenan_work(im,n_dims,sz_im,imout);

  return;
}


/*
 * Do the dilation. This works by first assembling a list of the pixels that have value NaN. Then, NaN-pixels with non-NaN nearest neighbors are replaced by the average of these neighbors.
 * One options would be to do this in a way that prioritizes pixels with the most neighbors; only those with the most neighbors are done (all in parallel). Afterwards, the number of nearest neighbors would be reassessed. However, let's try the simplest thing first.
 */
void imdilatenan_work(const float *im,int n_dims,const int *sz_im,float *imout)
{
  int *imneighbor_skip,*coords;
  long *nanlist,*nanlist_end,*nanlist_tmp,*nanlist_cur;
  float *newvalue,*newvalue_end;
  int dimIterator,n_pixels,neighbornum;
  long pixelIndex;
  float neighborsum,tmp;

  // Compute the pointer offset between adjacent pixels along all dimensions
  imneighbor_skip = (int*) mxMalloc(n_dims*sizeof(int));
  imneighbor_skip[0] = 1; 
  for (dimIterator = 1; dimIterator < n_dims; dimIterator++)
    imneighbor_skip[dimIterator] = imneighbor_skip[dimIterator-1] * sz_im[dimIterator-1];
  n_pixels = imneighbor_skip[n_dims-1] * sz_im[n_dims-1];
  
  // Allocate storage
  coords = (int*) mxMalloc(n_dims*sizeof(int));
  nanlist = (long*) mxMalloc(n_pixels * sizeof(long));
  nanlist_end = nanlist;
  newvalue = (float*) mxMalloc(n_pixels * sizeof(float));
  newvalue_end = newvalue;

  // Copy image data
  memcpy(imout,im,n_pixels*sizeof(float));

  // Create list of NaN pixels
  for (pixelIndex = 0; pixelIndex < n_pixels; pixelIndex++)
    if (mxIsNaN(imout[pixelIndex])) {
      *nanlist_end = pixelIndex;
      nanlist_end++;
    }

  // Go through list of NaNs
  while (nanlist_end != nanlist) {
    for (nanlist_tmp = nanlist, newvalue_end = newvalue; nanlist_tmp < nanlist_end; nanlist_tmp++, newvalue_end++) {
      // Calculate coordinates, so we know whether we're at an edge of
      // image (to decide which neighbors to examine)
      pixelIndex = *nanlist_tmp;
      for (dimIterator = n_dims-1; dimIterator >= 0; dimIterator--) {
	coords[dimIterator] = pixelIndex/imneighbor_skip[dimIterator];
	pixelIndex -= coords[dimIterator]*imneighbor_skip[dimIterator];
      }
      // Calculate the mean value of the non-NaN nearest neighbors
      // (along coordinate axes)
      neighborsum = 0;
      neighbornum = 0;
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
	if (coords[dimIterator] > 0)
	  if (!mxIsNaN(tmp = *(imout + *nanlist_tmp - imneighbor_skip[dimIterator]))) {
	    neighborsum += tmp;
	    neighbornum++;
	  }
	if (coords[dimIterator] < sz_im[dimIterator]-1)
	  if (!mxIsNaN(tmp = *(imout + *nanlist_tmp + imneighbor_skip[dimIterator]))) {
	    neighborsum += tmp;
	    neighbornum++;
	  }
      }
      if (neighbornum > 0)
	*newvalue_end = neighborsum/neighbornum; // average value
      else
	*newvalue_end = mxGetNaN();
    }
    // Replace values in the image
    for (nanlist_tmp = nanlist,newvalue_end = newvalue; nanlist_tmp < nanlist_end; nanlist_tmp++,newvalue_end++)
      imout[*nanlist_tmp] = *newvalue_end;
    // For the ones that were replaced with finite values, remove them
    // from the NaN list
    for (nanlist_tmp = nanlist,nanlist_cur = nanlist,newvalue_end=newvalue; nanlist_tmp < nanlist_end; nanlist_tmp++,newvalue_end++)
      if (mxIsNaN(*newvalue_end)) {
	// Keep it in the list
	*nanlist_cur = *nanlist_tmp;
	nanlist_cur++;
      }
    nanlist_end = nanlist_cur;
  }
  mxFree(coords);
  mxFree(nanlist);
  mxFree(newvalue);
}
