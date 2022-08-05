//#include <math.h>
//#include <string.h>
#include "mex.h"
#include "imiterators.cxx"
#include "image_utilities.h"

/* Syntax:
 *   imr = imrestrict(im)
 *   imr = imrestrict(im,dimFlag)
 */

// forward declarations
void imr_work(const float *im,const int n_dims,const int *sz_im,const bool* restrict_dim,float *imout,const int *sz_imout);

/*
 * This is the Matlab wrapper
 */

#define MAX_DIMS 3
#define THREE_TO_MAX_DIMS 27

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *im;
  float *imout;
  int n_dims,dimIndex;
  int sz_imout[MAX_DIMS];
  const int *sz_im;
  bool restrict_dim[MAX_DIMS];
  const mxArray *curarg;
  const bool* dimFlag;

  if (nrhs < 1 || nrhs > 2)
    mexErrMsgTxt("imrestrict: requires one or two inputs");
  if (nlhs != 1)
    mexErrMsgTxt("imrestrict: requires one output");

  // Parse the input
  // image
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("imrestrict: data must be a real single-precision array");
  im = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  if (n_dims > 3)
    mexErrMsgTxt("imrestrict: supports only images of dimensionality 3 or fewer");
  sz_im = mxGetDimensions(curarg);

  // dimFlag
  if (nrhs > 1) {
    curarg = prhs[1];
    if (!mxIsLogical(curarg))
      mexErrMsgTxt("imrestrict: dimFlag must be a real logical vector");
    if (mxGetNumberOfElements(curarg) != n_dims)
      mexErrMsgTxt("imrestrict: number of dimensions in image must match length of dimFlag");
    dimFlag = (bool*) mxGetData(curarg);
    // Copy it so can be modified (e.g., unity dimensions)
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      restrict_dim[dimIndex] = (dimFlag[dimIndex] > 0);
  } else
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      restrict_dim[dimIndex] = true;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (sz_im[dimIndex] == 1)
      restrict_dim[dimIndex] = false;



  // Calculate the size of the output
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (restrict_dim[dimIndex]) {
      sz_imout[dimIndex] = (sz_im[dimIndex]-1)/2;
      if (sz_imout[dimIndex] < 1)
	sz_imout[dimIndex] = 1;
    }
    else
      sz_imout[dimIndex] = sz_im[dimIndex];

  // Set up the output
  plhs[0] = mxCreateNumericArray(n_dims,sz_imout,
				 mxSINGLE_CLASS,mxREAL);
  imout = (float *) mxGetData(plhs[0]);

  // Do the actual work
  imr_work(im,n_dims,sz_im,restrict_dim,imout,sz_imout);

  return;
}

void imr_work(const float *im,const int n_dims,const int *sz_im,const bool* restrict_dim,float *imout,const int *sz_imout)
{
  float weights[THREE_TO_MAX_DIMS];
  int pixoffset[THREE_TO_MAX_DIMS];
  int coords[MAX_DIMS+1],coords_min[MAX_DIMS+1],coords_max[MAX_DIMS+1];
  long pixel_skip[MAX_DIMS];
  int dimIndex,nbrIndex,n_nbrs;
  float w_norm;
  const float *imptr;
  const int *pixoffset_ptr,*pixoffset_end;
  const float *weights_ptr;

  // Initialize the pixel offsets
  
  calc_pixel_skip(sz_im,n_dims,pixel_skip);
  /*
  mexPrintf("pixel_skip: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    mexPrintf("%d ",pixel_skip[dimIndex]);
  mexPrintf("\n");
  */

  // Calculate the weights and memory offsets for all of the nearest
  // neighbors; these are points for which the coordinate differences
  // are -1, 0, or 1. Dimensions that are not being restricted are
  // only allowed offsets of 0.
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    coords_min[dimIndex] = - int(restrict_dim[dimIndex]);
    coords_max[dimIndex] = int(restrict_dim[dimIndex] && sz_im[dimIndex] > 2); 
    coords[dimIndex] = coords_min[dimIndex];
  }
  // We use 1 extra coordinate as a sentinel to make
  // "carrying" style iteration easy.
  coords_min[n_dims] = -1;
  coords_max[n_dims] = 0;
  coords[n_dims] = -1;
  nbrIndex = 0;
  do {
    pixoffset[nbrIndex] = 0;
    weights[nbrIndex] = 1.0;
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
      pixoffset[nbrIndex] += coords[dimIndex] * pixel_skip[dimIndex];
      weights[nbrIndex] *= (1.0 + (coords[dimIndex] > coords_min[dimIndex] &&
				   coords[dimIndex] < coords_max[dimIndex]));
    }
    // Increment coords (with carry)
    dimIndex = 0;
    while (dimIndex <= n_dims && ++coords[dimIndex] > coords_max[dimIndex]) {
      coords[dimIndex] = coords_min[dimIndex];
      dimIndex++;
    }
    nbrIndex++;
  } while(coords[n_dims] < coords_max[n_dims]);
  n_nbrs = nbrIndex;
  pixoffset_end = pixoffset+n_nbrs;
  // Normalize the weights. (We calculate this from the data rather
  // than analytically to reduce the likelihood of getting confused by
  // the myriad possible cases.)
  w_norm = 0;
  for (nbrIndex = 0; nbrIndex < n_nbrs; nbrIndex++)
    w_norm += weights[nbrIndex];
  for (nbrIndex = 0; nbrIndex < n_nbrs; nbrIndex++)
    weights[nbrIndex] /= w_norm;
  /*
  for (nbrIndex = 0; nbrIndex < n_nbrs; nbrIndex++) {
    mexPrintf("pixoffset %d, weight %g\n",pixoffset[nbrIndex],weights[nbrIndex]);
  }
  */

  // Iterate over pixels of the output
  pixIterator pixI(sz_imout,n_dims,false);
  for (; !pixI.at_end(); pixI++,imout++) {
    /*
    mexPrintf("Coords: ");
    for (dimIndex = 0; dimIndex < pixI.nDims(); dimIndex++)
      mexPrintf("%d ",pixI.coord(dimIndex));
    mexPrintf("\n");
    */
    if (pixI.coord(0) == 0) {
      // Reset the pointer into the input image; the pointer is to the
      // 2j+1 index in each dimension
      imptr = im;
      for (dimIndex = 0; dimIndex < pixI.nDims(); dimIndex++)
	imptr += ((1+restrict_dim[dimIndex])*pixI.coord(dimIndex) +
		  restrict_dim[dimIndex]) * pixel_skip[dimIndex];
    }
    else
      imptr += (1+restrict_dim[0]);
    //mexPrintf("imptr %d\n",imptr-im);
    // Calculate the weighted sum of neighboring values
    *imout = 0;
    for (pixoffset_ptr = pixoffset,weights_ptr = weights;
	 pixoffset_ptr < pixoffset_end;
	 pixoffset_ptr++, weights_ptr++)
      *imout += imptr[*pixoffset_ptr] * *weights_ptr;
  }
}


