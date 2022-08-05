// Define MAIN while compiling to create a stand-alone command-line
// program. See Makefile.
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#endif

#include <vector>
//#include <math.h>
#include "image_utilities.h"
#include "timer_g.h"


//#define DEBUG_TIME 1

/*
 * Copyright 2010 by Timothy E. Holy
 */

/*
 * Implementation notes (see imqinterp.m for the call syntax):
 *
 * 1. A motivation for this syntax is to save the overhead of recalculating
 * the coefficients when multiple images need to be interpolated at the same
 * grid locations (e.g., the image and the weight mask).  So the
 * calculation of a single coefficient can be done based on g, and then there
 * should be a loop to apply that coefficient to the different images.
 * One should also think carefully about the bounds check, as this is the
 * most expensive part in the current matlab code.
 * Sketching it out: one can think of g as being n_pixels_out-by-n_dims_in
 * (even if this isn't the right shape), so all you really need to know is
 * the # of output pixels.  So the "work function" (the thing that does
 * the real work) is something like this:
 *    Prepare ahead: the memory offset list
 *    Loop over output pixels of g {  // good place for multithreading!
 *      Round & compute fractional components of g(x)
 *      Check bounds; if we are over the edge, (1) set this pixel to
 *        NaN in the outputs, and (2) go on to the next output pixel
 *      Compute the coefficients for all 3^n_dims_in terms for the
 *        interpolation, and all n_dims_in*3^n_dims_in for the
 *        gradients (if applicable), by looping over the neighbor
 *        combinations (-1, 0, 1). Store all of these in the
 *        coefficient list.
 *      Loop over input images and use these coefficients
 *    }
 *
 * 
 * 2. I would handle the variety in inputs by making a vector of pointers,
 * each corresponding to a single-channel image (so an RGB image would end up
 * getting represented as 3 pointers).  You can make a single list for all
 * the varargin inputs.  Similarly for the outputs.  That way the loop
 * over input images does not have to know anything about how these are
 * configured.
 *
 */ 

// Maximum number of dimensions is hardcoded to make life easy, could
// be done dynamically.
#define DIMS_MAX 3  // number of spatial dimensions
#define DIMS_MAX_TOTAL DIMS_MAX+1  // includes extra dimensions like RGB color
#define THREE_TO_DIMS_MAX 27   // 3^DIMS_MAX
#define OFFSET 1  // matlab is unit-offset

template <class T>
void mexWrapper(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[]);
template <class T>
void imqinterpwork(const T* g[],int n_dims_in,const int sz_out[],int n_images,const T** im,const int sz_im[],const int z[THREE_TO_DIMS_MAX][DIMS_MAX],const long im_memoffset[THREE_TO_DIMS_MAX],T** imout,T** grad_out,T badPixelValue);
int prepare_z(int z[THREE_TO_DIMS_MAX][DIMS_MAX],int n_dims);
template <class T>
long cumprod(T *sz,int n_dims);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  mxClassID imclass;

#ifdef DEBUG_TIME
  start_timer();
  printf("enter mexFunction @ %.4f\n", read_timer());
#endif

  if (nrhs < 2)
    mexErrMsgTxt("imqinterp: requires at least two inputs");
  if (nlhs < 1)
    mexErrMsgTxt("imqinterp: requires at least one output");
  if (nlhs != nrhs-1 && nlhs != 2*(nrhs-1))
    mexErrMsgTxt("imqinterp: the number of input and output images do not match");
  // Determine the data type, then call the appropriate templated function
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("imqinterp: g must be a real array");
  imclass = mxGetClassID(curarg);
  if (imclass == mxSINGLE_CLASS)
    mexWrapper<float>(nlhs,plhs,nrhs,prhs);
  else if (imclass == mxDOUBLE_CLASS)
    mexWrapper<double>(nlhs,plhs,nrhs,prhs);
  else
    mexErrMsgTxt("imqinterp: data type must be single or double");

#ifdef DEBUG_TIME
  printf("leave mexFunction @ %.4f\n", read_timer());
#endif
}

template <class T>
void mexWrapper(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  mxClassID imclass;
  const T *g[DIMS_MAX];
  const T *tmpP;
  int n_dims_g, n_dims_in,n_dims_out,n_dims_tmp;
  int n_pixels_in,n_pixels_out,n_images,n_nbrs;
  int dimIndex,imIndex,imIndex2,outputDimIndex,i;
  const mwSize *sz_tmp;
  int sz_g[DIMS_MAX+1];
  int sz_in[DIMS_MAX_TOTAL],sz_im_tmp[DIMS_MAX_TOTAL];
  int sz_out[DIMS_MAX_TOTAL+1];  // +1 for gradient
  long pixel_skip[DIMS_MAX];
  int z[THREE_TO_DIMS_MAX][DIMS_MAX];
  long memoffset[THREE_TO_DIMS_MAX];
  long tmpL;
  std::vector<const T*> imin;
  std::vector<T*> imout,imgradout;

  bool calculating_gradient = (nlhs == 2*(nrhs-1));

  // Parse the inputs
  // g
  curarg = prhs[0];
  n_dims_g = mxGetNumberOfDimensions(curarg);
  if (n_dims_g > DIMS_MAX + 1)
    mexErrMsgTxt("imqinterp: g has too many dimensions (fix by recompiling)");
  sz_tmp = mxGetDimensions(curarg);
  // Now get the dimensions of g again, skipping unity dimensions
  n_dims_g = skip_unity_dimensions(sz_tmp,n_dims_g,sz_g,DIMS_MAX+1);
  imclass = mxGetClassID(curarg);
  if (n_dims_g > 1) {
    n_dims_in = sz_g[n_dims_g-1];
    n_dims_out = n_dims_g-1;
  }
  else {
    n_dims_in = 1;
    n_dims_out = 1;
  }
  // Prepare the spatial parts of the output size
  for (dimIndex = 0; dimIndex < n_dims_out; dimIndex++)
    sz_out[dimIndex] = sz_g[dimIndex];
  if (n_dims_in < n_dims_out)
    mexWarnMsgTxt("imqinterp: g indicates fewer dimensions in images than in the output");
  n_pixels_out = cumprod(sz_out,n_dims_out);
  tmpP = (T*) mxGetData(curarg);
  for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
    g[dimIndex] = tmpP + dimIndex*n_pixels_out;
  
  // images
  for (imIndex = 0; imIndex < nrhs-1; imIndex++) {
    curarg = prhs[1+imIndex];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
      mexErrMsgTxt("imqinterp: image must be a real array");
    if (mxGetClassID(curarg) != imclass)
      mexErrMsgTxt("imqinterp: images must be of the same class as g");
    n_dims_tmp = mxGetNumberOfDimensions(curarg);
    if (n_dims_tmp > DIMS_MAX_TOTAL)
      mexErrMsgTxt("imqinterp: image dimensionality is too large (fix by recompiling)");
    sz_tmp = mxGetDimensions(curarg);
    n_dims_tmp = skip_unity_dimensions(sz_tmp,n_dims_tmp,sz_im_tmp,DIMS_MAX_TOTAL);
    if (imIndex == 0) {
      // Store the size of the first image
      for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
	sz_in[dimIndex] = sz_im_tmp[dimIndex];
      n_pixels_in = calc_pixel_skip(sz_in,n_dims_in,pixel_skip);
    } else {
      // Compare to size of the first image
      for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
	if (sz_im_tmp[dimIndex] != sz_in[dimIndex])
	  mexErrMsgTxt("imqinterp: all images must have the same spatial size");
    }
    // Append any extra dimensions of the input to the output size
    for (outputDimIndex = n_dims_out, dimIndex = n_dims_in; dimIndex < n_dims_tmp; dimIndex++,outputDimIndex++)
      sz_out[outputDimIndex] = sz_im_tmp[dimIndex];
    if (n_dims_tmp > n_dims_in)
      n_images = cumprod(sz_im_tmp+n_dims_in,n_dims_tmp-n_dims_in);
    else
      n_images = 1;
    // Push pointers to individual input images onto the heap
    tmpP = (T*) mxGetData(curarg);
    for (imIndex2 = 0; imIndex2 < n_images; imIndex2++, tmpP += n_pixels_in)
      imin.push_back(tmpP);
    // Allocate output storage
    plhs[imIndex] = mxCreateNumericArray(outputDimIndex,sz_out,imclass,mxREAL);
    // Push pointers to individual output images onto the heap
    tmpP = (T*) mxGetData(plhs[imIndex]);
    for (imIndex2 = 0; imIndex2 < n_images; imIndex2++, tmpP += n_pixels_out)
      imout.push_back((T*) tmpP);
    // Do the same steps for the gradient. The main issue is deciding
    // how to store the components of the gradient; here the choice is
    // [spatial_coordinates gradient_coordinate extra_coordinates]
    if (calculating_gradient) {
      // Create the output size vector
      sz_out[n_dims_out] = n_dims_in;
      for (outputDimIndex = n_dims_out, dimIndex = n_dims_in; dimIndex < n_dims_tmp; dimIndex++,outputDimIndex++)
	sz_out[outputDimIndex] = sz_im_tmp[dimIndex];
      // Allocate
      plhs[imIndex+nrhs-1] = mxCreateNumericArray(outputDimIndex+1,sz_out,
						  imclass,mxREAL);
      // Push pointers. Use just one pointer for all n_dims_in
      // components of the gradient; the work function will use the
      // appropriate pointer arithmetic.
      tmpP = (T*) mxGetData(plhs[imIndex+nrhs-1]);
      for (imIndex2 = 0; imIndex2 < n_images; imIndex2++, tmpP += n_pixels_out*n_dims_in)
	imgradout.push_back((T*) tmpP);
      
    }
  }

  // Prepare z and the memory offsets
  n_nbrs = prepare_z(z,n_dims_in);
  for (i = 0; i < n_nbrs; i++) {
    tmpL = 0;
    for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
      tmpL = tmpL + z[i][dimIndex]*pixel_skip[dimIndex];
    memoffset[i] = tmpL;
  }

  // Do the actual work
  if (calculating_gradient)
    imqinterpwork(g,n_dims_in,sz_out,(int) imin.size(),&imin[0],sz_in,z,memoffset,(T**) &imout[0],&imgradout[0],(T) mxGetNaN());
  else
    imqinterpwork(g,n_dims_in,sz_out,(int) imin.size(),&imin[0],sz_in,z,memoffset,(T**) &imout[0],(T**) NULL,(T) mxGetNaN());
  
  return;
}


/*
 * Do the interpolation & gradient calculation.
 * Note: keep all memory allocation out of this function, so it can be generic & fast.
 */
template <class T>
void imqinterpwork(const T* g[],int n_dims_in,const int sz_out[],int n_images,const T** im,const int sz_im[],const int z[THREE_TO_DIMS_MAX][DIMS_MAX],const long im_memoffset[THREE_TO_DIMS_MAX],T** imout,T** grad_out,T badPixelValue)
{
  int dimIndex,coefIndex,imIndex;
  bool have_coefficients,over_edge;
  T val,alpha,tmp,coef;
  long inputPixelIndex,outputPixelIndex,n_pixels_out;
  long pixel_skip[DIMS_MAX];
  T coords_frac[DIMS_MAX];
  int coords_int[DIMS_MAX];
  T interp_coef1[DIMS_MAX][3], interp_coef[THREE_TO_DIMS_MAX], interp_coef_grad[DIMS_MAX][THREE_TO_DIMS_MAX];
										
  bool calculating_gradient = (grad_out != NULL);
  int n_pts_interp = 1;
  for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
    n_pts_interp *= 3;
  calc_pixel_skip(sz_im,n_dims_in,pixel_skip);

  // Calculate the total number of output pixels
  n_pixels_out = 1;
  for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
    n_pixels_out *= sz_out[dimIndex];

  //
  // Loop over pixels of g (& the output)
  //
  // (start multithreading here)
  for (outputPixelIndex = 0; outputPixelIndex < n_pixels_out; outputPixelIndex++) {
    // Round & compute fractional components of g(x)
    for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++) {
      coords_int[dimIndex] = (int)(g[dimIndex][outputPixelIndex] + 0.5);
      coords_frac[dimIndex] = g[dimIndex][outputPixelIndex] - coords_int[dimIndex];
      coords_int[dimIndex] -= OFFSET;
    }
    // Test to see if this output pixel uses any input pixels that
    // are over-the-edge
    over_edge = false;
    for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
      if (coords_int[dimIndex] < 1 || 
	  coords_int[dimIndex] >= sz_im[dimIndex]-1)
	over_edge = true;
    if (over_edge) {
      for (imIndex = 0; imIndex < n_images; imIndex++)
	imout[imIndex][outputPixelIndex] = badPixelValue;
      continue;  // go on to the next output pixel
    }
    // Determine the memory offset of the "base" input point
    inputPixelIndex = 0;
    for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
      inputPixelIndex += coords_int[dimIndex]*pixel_skip[dimIndex];
    //
    // Loop over "output components": -1 = interpolated value,
    //   0 and higher = component of the gradient.
    //
    for (int outputComponent = -1; outputComponent < calculating_gradient*n_dims_in; outputComponent++) {
      //
      // Calculate interpolation coefficients
      //
      // Form the coefficients for each dimension
      for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++) {
	alpha = coords_frac[dimIndex];
	if (dimIndex == outputComponent) {
	  // this is a gradient component
	  interp_coef1[dimIndex][0] = alpha-0.5;
	  interp_coef1[dimIndex][2] = alpha+0.5;
	  interp_coef1[dimIndex][1] = -2*alpha;
	} else {
	  // this is either a value component, or a gradient component
	  // but not the dimension of the current derivative
	  tmp = alpha-0.5;
	  interp_coef1[dimIndex][0] = tmp*tmp/2;
	  tmp = alpha+0.5;
	  interp_coef1[dimIndex][2] = tmp*tmp/2;
	  interp_coef1[dimIndex][1] = 0.75 - alpha*alpha;
	}
      }
      // Compute products of the one-dimensional coefficients
      for (coefIndex = 0; coefIndex < n_pts_interp; coefIndex++) {
	for (coef = 1, dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
	  coef *= interp_coef1[dimIndex][1+z[coefIndex][dimIndex]];
	interp_coef[coefIndex] = coef;
      }
      //
      // Loop over the different images
      //
      for (imIndex = 0; imIndex < n_images; imIndex++) {
	// Calculate the interpolated value
	val = 0;
	for (coefIndex = 0; coefIndex < n_pts_interp; coefIndex++)
	  val += interp_coef[coefIndex] * 
	    im[imIndex][inputPixelIndex+im_memoffset[coefIndex]];
	if (outputComponent < 0)
	  imout[imIndex][outputPixelIndex] = val;  // interpolated image
	else
	  grad_out[imIndex][outputPixelIndex+outputComponent*n_pixels_out] = val;  // gradient component
      }  // loop over images
    }  // loop over outputComponent
  }  // loop over outputPixel
}


/* 
 * A utility function to prepare all coordinate offsets from the
 * "current" pixel, using the -1, 0, and +1 neighbors along each axis.
 */
int prepare_z(int z[][DIMS_MAX],int n_dims)
{
  int dimIndex,pixelIndex,n_nbrs;

  n_nbrs = 1;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    z[0][dimIndex] = -1;
    n_nbrs *= 3;  // becomes 3^n_dims
  }

  for (pixelIndex = 1; pixelIndex < n_nbrs; pixelIndex++) {
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      z[pixelIndex][dimIndex] = z[pixelIndex-1][dimIndex];
    // Implement add one with carry
    dimIndex = 0;
    z[pixelIndex][0] = z[pixelIndex-1][0]+1;
    while (z[pixelIndex][dimIndex] > 1 && dimIndex < n_dims-1) {
      z[pixelIndex][dimIndex] = -1;
      dimIndex++;
      z[pixelIndex][dimIndex] = z[pixelIndex-1][dimIndex]+1;  // the carry
    }
  }
  return n_nbrs;
}

template <class T>
long cumprod(T *sz,int n_dims)
{
  long n_pixels_out = 1;
  for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
    n_pixels_out *= sz[dimIndex];
  return n_pixels_out;
}



#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 2;
  const int n_outputs = 2;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "g",
    "im"
  };
  const char *output_names[] = {
    "imo",
    "grado"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables g and im\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
