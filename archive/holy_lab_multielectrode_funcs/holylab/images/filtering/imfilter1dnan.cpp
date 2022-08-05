#include "mex.h"
#include <math.h>

//#if defined(WIN) || defined(_WIN)
//#include <float.h>
//#define isnan _isnan
//#endif


/*
 * imfilter1dnan: filter an image along one dimension, ignoring fill values
 *
 * Syntax:
 *   imout = imfilter1dnan(im,h,dim)
 *   imout = imfilter1dnan(im,h,dim,stepsize)
 * where
 *   im is the input image (may be multidimensional)
 *   h is the 1-d filter
 *   dim is the dimension to filter along
 * and
 *   imout is the filtered image.
 *
 * Copyright 2006 by Timothy E. Holy
 */

void imfilter1dnan_work(const float *im,int n_dims,const int *sz_im,int *coords,const double *h,int lh,int dimIndex,int stepsize,float *imout);
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
  const double *h;
  float *imout;//,fillval;
  int n_dims,dimIndex,lh,stepsize;
  const int *sz_im,*sz_h;
  int *coords;
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("imfilter1dnan: requires three or four inputs");
  if (nlhs != 1)
    mexErrMsgTxt("imfilter1dnan: requires one output");

  // Parse the inputs
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("imfilter1dnan: image must be a real single-precision array");
  im = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz_im = mxGetDimensions(curarg);

  // filter (h)
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("imfilter1dnan: filter must be a real double-precision array");
  h = mxGetPr(curarg);
  if (!isVector(mxGetDimensions(curarg),mxGetNumberOfDimensions(curarg)))
    mexErrMsgTxt("imfilter1dnan: filter (h) must be a vector");
  lh = mxGetNumberOfElements(curarg);
  if (lh != 2*((lh-1)/2)+1)
    mexErrMsgTxt("imfilter1dnan: filter must have odd length");

  // dimIndex
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("imfilter1dnan: dim must be a real (double) integer");
  if (mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("imfilter1dnan: dim must be a scalar");
  dimIndex = (int) mxGetScalar(curarg) - 1;  // switch to zero offset
  if (dimIndex < 0 || dimIndex >= n_dims)
    mexErrMsgTxt("imfilter1dnan: dim must be a valid dimension index");


  // stepsize
  stepsize = 1;
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
      mexErrMsgTxt("imfilter1dnan: stepsize must be a real numeric value");
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("imfilter1dnan: stepsize must be a scalar");
    stepsize = (int) mxGetScalar(curarg);
    if (stepsize < 1)
      mexErrMsgTxt("imfilter1dnan: stepsize must be at least 1");
  }

  // Create working storage
  coords = (int*) mxMalloc(n_dims*sizeof(int));

  // Set up the output
  plhs[0] = mxCreateNumericArray(n_dims,sz_im,mxSINGLE_CLASS,mxREAL);
  imout = (float *) mxGetData(plhs[0]);

  // Do the actual work
  imfilter1dnan_work(im,n_dims,sz_im,coords,h,lh,dimIndex,stepsize,imout);

  mxFree(coords);

  return;
}


/*
 * Do the filtering.
 */
void imfilter1dnan_work(const float *im,int n_dims,const int *sz_im,int *coords,const double *h,int lh,int dimIndex,int stepsize,float *imout)
{
  const float *neighborpixel;
  double pixel_sum,h_sum;
  int filter_halfsize,imneighbor_skip,n_pixels,n_remaining;
  int dimIterator,pixelIterator,filterIterator,filterStart,filterStop;
  int continuing;

  // Shift the filter pointer to the center of the filter, so we
  // index in both + and - directions
  filter_halfsize = (lh-1)/2;
  h += filter_halfsize;

  // Compute the pointer offset between adjacent pixels along the chosen dimension
  imneighbor_skip = 1;
  dimIterator = 0;
  while (dimIterator < dimIndex)
    imneighbor_skip *= sz_im[dimIterator++];
  // Continue onward to compute the total number of pixels in the image
  n_pixels = imneighbor_skip;
  while (dimIterator < n_dims)
    n_pixels *= sz_im[dimIterator++];
  // Adjust imneighbor_skip to accomodate the possibility that we'll be skipping over neighbors
  imneighbor_skip *= stepsize;

  // Initialize coords
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
    coords[dimIterator] = 0;
  //mexPrintf("Fill value: %g\n",fillval);
  
  // Iterate over pixels in the image. This is reasonably cache-efficient, because when we move on to the next pixel, the cache lines should typically be ready for the next pixel.
  for (pixelIterator = 0; pixelIterator < n_pixels; pixelIterator++, im++, imout++) {
    //mexPrintf("%d: %f, isnan: %d \n",pixelIterator,*im,isnan(*im));
    if (mxIsNaN(*im)) {
      *imout = *im;
      //mexPrintf("%d: Filling\n",pixelIterator);
    }
    else {
      pixel_sum = h_sum = 0;
      // Determine the subset of the filter that is within the boundaries of the image
      filterStart = -filter_halfsize;
      if (coords[dimIndex] - stepsize * filter_halfsize < 0)
	filterStart = (int) ceil(-(double)coords[dimIndex]/stepsize);
      filterStop = filter_halfsize;
      if (coords[dimIndex] + stepsize * filter_halfsize >= sz_im[dimIndex])
	filterStop = (int) floor(((double)sz_im[dimIndex] - 1.0 - coords[dimIndex])/stepsize);
      //n_remaining = sz_im[dimIndex] - coords[dimIndex] - 1;
      //if (n_remaining < filterStop)
      //filterStop = n_remaining;
      //mexPrintf("%d: start %d, stop %d\n",pixelIterator,filterStart,filterStop);
      // Filter the image pixels
      for (filterIterator = filterStart, neighborpixel = im + filterStart*imneighbor_skip; filterIterator <= filterStop; filterIterator++, neighborpixel += imneighbor_skip) {
	//mexPrintf("%g ",*neighborpixel);
	//if (*neighborpixel != fillval) && !isnan(*neighborpixel)) {
	if (!mxIsNaN(*neighborpixel)) {
	  pixel_sum += h[filterIterator] * (*neighborpixel);
	  h_sum += h[filterIterator];
	}
      }
      //mexPrintf("\n");
      // Save the result
      *imout = pixel_sum/h_sum;
    }
    // Update the coordinate counters
    dimIterator = 0;
    while (dimIterator < n_dims) {
      if (++coords[dimIterator] >= sz_im[dimIterator]) {
	coords[dimIterator] = 0;
	dimIterator++;
      }
      else
	break;
    }
  }
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
