#include "mex.h"
#include <math.h>
#include "image_utilities.h"

/*
 * iminterp: interpolate an image quickly
 *
 * Syntax:
 *   imout = iminterp(im,X1,X2,...)
 *   [imout,w] = iminterp(im,X1,X2,...)
 *   imout = iminterp(im,X1,X2,...,'extrap')
 * where
 *   im is the input image (may be multidimensional)
 *   X1, X2, ... are arrays of the same dimensions as the image,
 *     specifying a set of evaluation locations.  These arrays specify a
 *     deformation X1 = g_1(x) for x over the points on a grid. Note:
 *     any dimensions of size 1 (e.g., for a 1-by-100 "image") are
 *     ignored, and should not have corresponding Xi.
 * and
 *   imout is the interpolated image
 *   w is a weight array with the same size as the image.  w is 1 in
 *     the interior but between 0 and 1 on the boundary.  For example,
 *     in d=1, a pixel evaluated at x = 0.8 is just beyond the left edge
 *     of the image.  Without the w output, the interpolated value would
 *     be NaN.  With the w output, the pixel would have the value of the
 *     pixel at x=1 but would have a weight of 0.8.  This is useful in
 *     ensuring continuity of functions of interpolated images (e.g., in
 *     image registration).
 *
 * If you want to have the image value linearly extrapolated beyond
 * the edge of the image, supply the extra string argument
 * 'extrap'. This mode is incompatible with the "w" weight array
 * output.
 *
 * Copyright 2006-2009 by Timothy E. Holy
 *  2009: support for double added (TEH)
 */

#define DIMS_MAX 3

template <class T>
void iminterpwork(const T *im,int n_dims_im,const int *sz_im,const void *g[],int n_dims_g,const int *sz_g,T *imout,T *w,const bool extrapflag);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const void *im,*g[DIMS_MAX];
  void *imout,*w;
  int n_dims,n_dims_im,n_dims_g,n_dims_tmp,dimIndex;
  const int *sz_im_tmp,*sz_g_tmp;
  int sz_im[DIMS_MAX], sz_g[DIMS_MAX], sz_g_tst[DIMS_MAX];
  const mxArray *curarg;
  bool extrapflag;
  mxClassID imclass;

  if (nrhs < 2)
    mexErrMsgTxt("iminterp: requires at least two inputs");
  if (nlhs < 1 || nlhs > 2)
    mexErrMsgTxt("iminterp: requires one or two outputs");

  // Parse the inputs
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("iminterp: image must be a real array");
  imclass = mxGetClassID(curarg);
  n_dims_tmp = mxGetNumberOfDimensions(curarg);
  sz_im_tmp = mxGetDimensions(curarg);
  // Don't use Matlab's notion of dimensionality, because otherwise we
  // have problems with 1-dimensional "images"
  n_dims_im = skip_unity_dimensions(sz_im_tmp,n_dims_tmp,sz_im,DIMS_MAX);
  if (n_dims_im < 1)
    mexErrMsgTxt("iminterp: supports only dimensions 1, 2, or 3");
  im = mxGetData(curarg);
  if (n_dims_im != nrhs-1 && n_dims_im != nrhs-2)
    mexErrMsgTxt("iminterp: dimensionality of image and the number of coordinates do not match");

  // Components of g (X1, X2, ...)
  n_dims = n_dims_im;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    curarg = prhs[1+dimIndex];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
      mexErrMsgTxt("iminterp: components of g must be a real array");
    if (mxGetClassID(curarg) != imclass)
      mexErrMsgTxt("iminterp: class of g must match that of the image");
    g[dimIndex] = mxGetData(curarg);
    n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				       mxGetNumberOfDimensions(curarg),
				       sz_g_tst,DIMS_MAX);
    if (dimIndex == 0) {
      // Copy over g size info
      n_dims_g = n_dims_tmp;
      for (int dimIndex2 = 0; dimIndex2 < n_dims_g; dimIndex2++)
	sz_g[dimIndex2] = sz_g_tst[dimIndex2];
    }
    else
      if (!validate_dimensions(sz_g_tst,n_dims_tmp,sz_g,n_dims_g))
	mexErrMsgTxt("iminterp: all components of g (X) must have the same size");
  }
  // Note: we don't enforce that dim(X1) = dim(im), because one might
  // want to use this to create a "slice" from an image

  extrapflag = false;
  if (n_dims == nrhs-2) {
    curarg = prhs[nrhs-1];
    if (!mxIsChar(curarg))
      mexErrMsgTxt("iminterp: extrapolation flag must be a string");
    if (!mxIsEmpty(curarg))
      if (mxGetScalar(curarg) == double('e')) {
	extrapflag = true;
      }
  }
  if (extrapflag && nlhs > 1)
    mexErrMsgTxt("iminterp: extrapolation is inconsistent with weight output");

  // Set up the output, an array of the size of X1
  plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),
				 mxGetDimensions(prhs[1]),
				 imclass,mxREAL);
  imout = mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),
				   mxGetDimensions(prhs[1]),
				   imclass,mxREAL);
    w = mxGetData(plhs[1]);
  }
  else
    w = NULL;

  // Do the actual work
  if (imclass == mxSINGLE_CLASS)
    iminterpwork<float>((float *)im,n_dims_im,sz_im,g,n_dims_g,sz_g,(float *)imout,(float *)w,extrapflag);
  else if (imclass == mxDOUBLE_CLASS)
    iminterpwork<double>((double *)im,n_dims_im,sz_im,g,n_dims_g,sz_g,(double *)imout,(double *)w,extrapflag);
  else
    mexErrMsgTxt("iminterp: image class not yet supported.");

  return;
}


/*
 * Do the interpolation.
 */
template <class T>
void iminterpwork(const T *im,int n_dims_im,const int *sz_im,const void *g_in[],int n_dims_g,const int *sz_g,T *imout,T *w,const bool extrapflag)
{
  long n_pixels_out,pixel_step[DIMS_MAX],pixel_offset;
  int dimIndex,coords_int[DIMS_MAX],use_w;
  T *imout_end;
  const T *current_pixel_ptr;
  T coords_frac[DIMS_MAX],weight,this_coord,this_deriv;
  bool is_outside,within_boundary;
  const T *g[DIMS_MAX];

  use_w = (w != NULL);

  // Cast g appropriately
  for (dimIndex = 0; dimIndex < n_dims_g; dimIndex++)
    g[dimIndex] = (T*) g_in[dimIndex];

  // Figure out when to stop
  n_pixels_out = 1;
  for (dimIndex = 0; dimIndex < n_dims_g; dimIndex++)
    n_pixels_out *= sz_g[dimIndex];
  imout_end = imout+n_pixels_out;

  // Determine pixel spacing in input image
  calc_pixel_skip(sz_im,n_dims_im,pixel_step);

  // Loop over pixels of the output
  // (start multithreading here)
  for (; imout < imout_end; imout++, w+=use_w) {
    // Calculate the pointer to the "floor" point, and the fraction
    // along each between-integer coordinate in each dimension
    pixel_offset = 0;
    is_outside = false;
    within_boundary = true;
    for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++) {
      coords_int[dimIndex] = int(floorf(*(g[dimIndex])));
      coords_frac[dimIndex] = *(g[dimIndex]) - coords_int[dimIndex];
      if (coords_int[dimIndex] < 1 || 
	  coords_int[dimIndex] >= sz_im[dimIndex])
	is_outside = true;
      if (is_outside && use_w)
	if (coords_int[dimIndex] < 0 ||
	    coords_int[dimIndex] >= sz_im[dimIndex]+1)
	  within_boundary = false;
    }
    //mexPrintf("%d %d is_outside: %d within_boundary: %d\n",coords_int[0],coords_int[1],is_outside,within_boundary);
    if (is_outside)
      if (extrapflag || (use_w && within_boundary)) {
	// Point is on the boundary. Use the nearest-neighbor
	// point for the image value, but adjust the weight to encode
	// the degree to which it's sticking beyond the data.
	// Or, extrapolate its value.
	weight = 1.0;
	for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++) {
	  this_coord = *(g[dimIndex]);
	  if (this_coord > 1 && this_coord < sz_im[dimIndex])
	    coords_int[dimIndex] = int(roundf(this_coord));
	  else if (this_coord <= 1) {
	    // this_coord is between 0 and 1, i.e., on the "left" boundary
	    coords_int[dimIndex] = 1;
	    weight *= this_coord;  // = coords_frac[dimIndex]
	  }
	  else {
	    // this_coord is between sz_im[dimIndex] and sz_im[dimIndex]+1,
	    // i.e., on the "right" boundary
	    coords_int[dimIndex] = sz_im[dimIndex];
	    weight *= 1.0 - coords_frac[dimIndex];
	  }
	  //mexPrintf("coords_int[%d] = %d\n",dimIndex,coords_int[dimIndex]);
	}
	for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++)
	  pixel_offset += (coords_int[dimIndex]-1)*pixel_step[dimIndex];
	*imout = im[pixel_offset];
	if (use_w)
	  *w = weight;
	if (extrapflag) {
	  // Correct the value by linear extrapolation
	  for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++) {
	    this_coord = *(g[dimIndex]) - coords_int[dimIndex];
	    if (coords_int[dimIndex] == 1)
	      this_deriv = im[pixel_offset+pixel_step[dimIndex]] - im[pixel_offset];
	    else if (coords_int[dimIndex] == sz_im[dimIndex])
	      this_deriv = im[pixel_offset] - im[pixel_offset-pixel_step[dimIndex]];
	    else
	      this_deriv = (im[pixel_offset+pixel_step[dimIndex]] -
			    im[pixel_offset-pixel_step[dimIndex]])/2;
	    /*
	    mexPrintf("coords_int %d, sz_im %d\n",coords_int[dimIndex],sz_im[dimIndex]);
	    mexPrintf("this_coord %g, this_deriv %g\n",this_coord,this_deriv);
	    */
	    *imout += this_coord * this_deriv;
	  }
	}
      }
      else {
	// Point was outside the boundary, or it's not in the interior
	// and we're not using w
	*imout = mxGetNaN();
	if (use_w)
	  *w = 0;
      }
    else {
      // Point is in the interior, we do the "regular" calculation
      if (use_w)
	*w = 1;
      for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++)
	pixel_offset += (coords_int[dimIndex]-1)*pixel_step[dimIndex];
      current_pixel_ptr = im + pixel_offset;
      // Do the linear interpolation
      if (n_dims_im == 1) {   // d = 1
	*imout = (1-coords_frac[0]) * (*current_pixel_ptr) +
	  coords_frac[0] * *(current_pixel_ptr+pixel_step[0]);
      } else if (n_dims_im == 2) {   // d = 2
	*imout = (1-coords_frac[0]) * (1-coords_frac[1]) * (*current_pixel_ptr) +
	  coords_frac[0] * (1-coords_frac[1]) * *(current_pixel_ptr+pixel_step[0]) +
	  (1-coords_frac[0]) * coords_frac[1] * *(current_pixel_ptr+pixel_step[1]) +
	  coords_frac[0] * coords_frac[1] * *(current_pixel_ptr+pixel_step[0]+pixel_step[1]);
      } else {              // d = 3
	*imout = (1-coords_frac[0]) * (1-coords_frac[1]) * (1-coords_frac[2]) * (*current_pixel_ptr) +
	  coords_frac[0] * (1-coords_frac[1]) * (1-coords_frac[2]) * *(current_pixel_ptr+pixel_step[0]) +
	  (1-coords_frac[0]) * coords_frac[1] * (1-coords_frac[2]) * *(current_pixel_ptr+pixel_step[1]) +
	  (1-coords_frac[0]) * (1-coords_frac[1]) * coords_frac[2] * *(current_pixel_ptr+pixel_step[2]) +
	  coords_frac[0] * coords_frac[1] * (1-coords_frac[2]) * *(current_pixel_ptr+pixel_step[0]+pixel_step[1]) +
	  (1-coords_frac[0]) * coords_frac[1] * coords_frac[2] * *(current_pixel_ptr+pixel_step[1]+pixel_step[2]) +
	  coords_frac[0] * (1-coords_frac[1]) * coords_frac[2] * *(current_pixel_ptr+pixel_step[0]+pixel_step[2]) +
	  coords_frac[0] * coords_frac[1] * coords_frac[2] * *(current_pixel_ptr+pixel_step[0]+pixel_step[1]+pixel_step[2]);
      }
    }
    // Move the "lookup table" forward in preparation for the next pixel
    for (dimIndex = 0; dimIndex < n_dims_im; dimIndex++)
      (g[dimIndex])++;
  }
  
}

