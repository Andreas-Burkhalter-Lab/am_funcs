// Define MAIN while compiling to create a stand-alone command-line
// program. See Makefile.
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
//#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#endif

#include <vector>
#include <unistd.h>
#include "image_utilities.h"
#include "threadpool.hpp"
#include "timer_g.h"

//#define DEBUG_TIME 1

/*
 * Copyright 2010 by Timothy E. Holy & Zhongsheng Guo
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

///thread related:
const int cMaxThreadNums=24; //the max number of threads
int nRunningThreads=0; //must <=cMaxThreadNums
ThreadPool pool(cMaxThreadNums);


template <class T>
void mexWrapper(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[]);
template <class T>
void imqinterpwork(const T* g[],int n_dims_in,const int sz_out[],int n_dims_out,int n_images,const T** im,const int sz_im[],const int z[THREE_TO_DIMS_MAX][DIMS_MAX],const long im_memoffset[THREE_TO_DIMS_MAX],T** imout,T** grad_out,T badPixelValue);
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
  printf("enter mexFuntion() @  %.4f\n", read_timer());
#endif

  if(!nRunningThreads){
#ifdef DEBUG_TIME
     printf("auto set nRunningThreads\n");
#endif
     nRunningThreads=sysconf(_SC_NPROCESSORS_ONLN)/2;
     if(nRunningThreads<1)nRunningThreads=1;
     if(nRunningThreads>4)nRunningThreads=4;
  }

  if(nrhs==1){
     /* The input must be a noncomplex scalar double.*/
     int mrows = mxGetM(prhs[0]);
     int ncols = mxGetN(prhs[0]);
     if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
         !(mrows==1 && ncols==1) ) {
        mexErrMsgTxt("Input must be a non-complex scalar double.");
     }

     // get #threads:
     nRunningThreads = (int) *mxGetPr(prhs[0]);
     if(nRunningThreads>cMaxThreadNums)nRunningThreads=cMaxThreadNums;
     return;
  }//if, setting #running threads only

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
  printf("leave mexFuntion() @  %.4f\n", read_timer());
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
  // Determining the dimensionality is a bit subtle. Difficult cases:
  //  one-dimensional "images", e.g., an image of size [5 1]. 
  //  asking for a scalar output from an image (e.g., g is of size [1 1 2])
  n_dims_g = mxGetNumberOfDimensions(curarg);
  if (n_dims_g > DIMS_MAX + 1)
    mexErrMsgTxt("imqinterp: g has too many dimensions (fix by recompiling)");
  sz_tmp = mxGetDimensions(curarg);
  n_dims_in = sz_tmp[n_dims_g-1];
  // Now get the dimensions of g again, skipping unity dimensions
  n_dims_out = skip_unity_dimensions(sz_tmp,n_dims_g-1,sz_g,DIMS_MAX+1);
  //mexPrintf("n_dims_g %d, n_dims_in %d, n_dims_out %d\n",n_dims_g,n_dims_in,n_dims_out);
  if (n_dims_out == 0) {
    n_dims_out = 1;  // user just wants value at a single point
    sz_g[0] = 1;
  }
  imclass = mxGetClassID(curarg);
  //if (n_dims_g > 0) {
  //  n_dims_in = sz_g[n_dims_g-1];
  //  n_dims_out = n_dims_g-1;
  //}
  //else {
  //  n_dims_in = 1;
  //  n_dims_out = 1;
  //}

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
    imqinterpwork(g,n_dims_in,sz_out,n_dims_out,(int) imin.size(),&imin[0],sz_in,z,memoffset,(T**) &imout[0],&imgradout[0],(T) mxGetNaN());
  else
    imqinterpwork(g,n_dims_in,sz_out,n_dims_out,(int) imin.size(),&imin[0],sz_in,z,memoffset,(T**) &imout[0],(T**) NULL,(T) mxGetNaN());
}


//NOTE: pxl range [outputPixelIndexFrom, outputPixelIndexTo)
template<typename T>
void procPxlRange(long outputPixelIndexFrom, long outputPixelIndexTo, int n_dims_in, const int sz_im[],
                  int n_images, const long pixel_skip[], bool calculating_gradient, int n_pts_interp, const int z[][DIMS_MAX],
                  const long im_memoffset[], long n_pixels_out,
                  const T* g[], T badPixelValue, T** imout,
                  const T** im, T** grad_out) {
   int coords_int[DIMS_MAX];
   T coords_frac[DIMS_MAX];
   T interp_coef1[DIMS_MAX][3], interp_coef[THREE_TO_DIMS_MAX]; //, interp_coef_grad[DIMS_MAX][THREE_TO_DIMS_MAX];


   for(long outputPixelIndex = outputPixelIndexFrom; outputPixelIndex < outputPixelIndexTo; outputPixelIndex++) {
      // Round & compute fractional components of g(x)
      for(int dimIndex = 0; dimIndex < n_dims_in; dimIndex++) {
         coords_int[dimIndex] = (int) (g[dimIndex][outputPixelIndex] + 0.5);
         coords_frac[dimIndex] = g[dimIndex][outputPixelIndex]
               - coords_int[dimIndex];
         coords_int[dimIndex] -= OFFSET;
      }
      // Test to see if this output pixel uses any input pixels that
      // are over-the-edge
      bool over_edge = false;
      for(int dimIndex = 0; dimIndex < n_dims_in; dimIndex++){
         if (coords_int[dimIndex] < 1 || coords_int[dimIndex] >= sz_im[dimIndex] - 1) over_edge = true;
      }
      if(over_edge) {
         for(int imIndex = 0; imIndex < n_images; imIndex++) imout[imIndex][outputPixelIndex] = badPixelValue;
         continue; // go on to the next output pixel
      }
      // Determine the memory offset of the "base" input point
      long inputPixelIndex = 0;
      for(int dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
         inputPixelIndex += coords_int[dimIndex] * pixel_skip[dimIndex];
      //
      // Loop over "output components": -1 = interpolated value,
      //   0 and higher = component of the gradient.
      //
      for(int outputComponent = -1; outputComponent < calculating_gradient * n_dims_in; outputComponent++) {
         //
         // Calculate interpolation coefficients
         //
         // Form the coefficients for each dimension
         for (int dimIndex = 0; dimIndex < n_dims_in; dimIndex++) {
            T alpha = coords_frac[dimIndex];
            if (dimIndex == outputComponent) {
               // this is a gradient component
               interp_coef1[dimIndex][0] = alpha - 0.5;
               interp_coef1[dimIndex][2] = alpha + 0.5;
               interp_coef1[dimIndex][1] = -2 * alpha;
            } else {
               // this is either a value component, or a gradient component
               // but not the dimension of the current derivative
               T tmp = alpha - 0.5;
               interp_coef1[dimIndex][0] = tmp * tmp / 2;
               tmp = alpha + 0.5;
               interp_coef1[dimIndex][2] = tmp * tmp / 2;
               interp_coef1[dimIndex][1] = 0.75 - alpha * alpha;
            }
         }
         // Compute products of the one-dimensional coefficients
         for (int coefIndex = 0; coefIndex < n_pts_interp; coefIndex++) {
            T coef = 1;
            for (int dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
               coef *= interp_coef1[dimIndex][1 + z[coefIndex][dimIndex]];
            interp_coef[coefIndex] = coef;
         }
         //
         // Loop over the different images
         //
         for (int imIndex = 0; imIndex < n_images; imIndex++) {
            // Calculate the interpolated value
            T val = 0;
            for (int coefIndex = 0; coefIndex < n_pts_interp; coefIndex++)
               val += interp_coef[coefIndex] * im[imIndex][inputPixelIndex + im_memoffset[coefIndex]];
            if (outputComponent < 0)  imout[imIndex][outputPixelIndex] = val; // interpolated image
            else  grad_out[imIndex][outputPixelIndex + outputComponent  * n_pixels_out] = val; // gradient component
         } // loop over images
      } // loop over outputComponent
   } // loop over outputPixel

}//procPxlRange(),

template <class T>
struct Params{
   long outputPixelIndexFrom;
   long outputPixelIndexTo;
   static int n_dims_in;
   static int* coords_int;
   static const int* sz_im;
   static int n_images;
   static long* pixel_skip;
   static bool calculating_gradient;
   static int n_pts_interp;
   static const int (*z)[DIMS_MAX];
   static const long *im_memoffset;
   static long n_pixels_out;
   static const T** g;
   static T* coords_frac;
   static T badPixelValue;
   static T** imout;
   static T (*interp_coef1)[3];
   static T* interp_coef;
   static const T** im;
   static T** grad_out;

   void setParams(long outputPixelIndexFrom, long outputPixelIndexTo){
      this->outputPixelIndexFrom=outputPixelIndexFrom;
      this->outputPixelIndexTo=outputPixelIndexTo;

   }

   void setCommonParams(int n_dims_in, const int sz_im[],
                  int n_images, long pixel_skip[], bool calculating_gradient, int n_pts_interp,  const int z[][DIMS_MAX],
                  const long im_memoffset[], long n_pixels_out,
                  const T* g[], T badPixelValue, T** imout,
                  const T** im, T** grad_out) {
      Params::n_dims_in=n_dims_in;
      Params::sz_im=sz_im;
      this->n_images=n_images;
      this->pixel_skip=pixel_skip;
      this->calculating_gradient=calculating_gradient;
      this->n_pts_interp=n_pts_interp;
      this->z=z;
      this->im_memoffset=im_memoffset;
      this->n_pixels_out=n_pixels_out;
      this->g=g;
      this->badPixelValue=badPixelValue;
      this->imout=imout;
      this->im=im;
      this->grad_out=grad_out;
   }//setParams(),
};//struct, Params
template<class T> int Params<T>::n_dims_in;
template<class T> int* Params<T>::coords_int;
template<class T> const int* Params<T>::sz_im;
template<class T> int Params<T>::n_images;
template<class T> long* Params<T>::pixel_skip;
template<class T> bool Params<T>::calculating_gradient;
template<class T> int Params<T>::n_pts_interp;
template<class T> const int (*Params<T>::z)[DIMS_MAX];
template<class T> const long *Params<T>::im_memoffset;
template<class T> long Params<T>::n_pixels_out;
template<class T> const T** Params<T>::g;
template<class T> T* Params<T>::coords_frac;
template<class T> T Params<T>::badPixelValue;
template<class T> T** Params<T>::imout;
template<class T> T (*Params<T>::interp_coef1)[3];
template<class T> T* Params<T>::interp_coef;
template<class T> const T** Params<T>::im;
template<class T> T** Params<T>::grad_out;

template <class T>
void threadEntry(void* arg){
   Params<T> * p=(Params<T> *)arg;
   procPxlRange(p->outputPixelIndexFrom, p->outputPixelIndexTo, p->n_dims_in, p->sz_im,
                p->n_images, p->pixel_skip, p->calculating_gradient, p->n_pts_interp, p->z,
                p->im_memoffset, p->n_pixels_out,
                p->g, p->badPixelValue, p->imout,
                p->im, p->grad_out);

}//threadEntry(),

/*
 * Do the interpolation & gradient calculation.
 * Note: keep all memory allocation out of this function, so it can be generic & fast.
 */
template <class T>
void imqinterpwork(const T* g[],int n_dims_in,const int sz_out[],int n_dims_out,int n_images,const T** im,const int sz_im[],const int z[THREE_TO_DIMS_MAX][DIMS_MAX],const long im_memoffset[THREE_TO_DIMS_MAX],T** imout,T** grad_out,T badPixelValue)
{
  int dimIndex,coefIndex,imIndex;
  bool have_coefficients;
  //bool over_edge;
  T val,alpha,tmp,coef;
  long inputPixelIndex,outputPixelIndex,n_pixels_out;
  long pixel_skip[DIMS_MAX];

  Params<T> params[cMaxThreadNums];
										
  bool calculating_gradient = (grad_out != NULL);
  int n_pts_interp = 1;
  for (dimIndex = 0; dimIndex < n_dims_in; dimIndex++)
    n_pts_interp *= 3;
  calc_pixel_skip(sz_im,n_dims_in,pixel_skip);

  // Calculate the total number of output pixels
  n_pixels_out = 1;
  for (dimIndex = 0; dimIndex < n_dims_out; dimIndex++)
    n_pixels_out *= sz_out[dimIndex];

  params[0].setCommonParams(n_dims_in, sz_im,
            n_images, pixel_skip, calculating_gradient, n_pts_interp, z,
            im_memoffset, n_pixels_out,
            g, badPixelValue, imout,
            im, grad_out);

   bool force_st=n_pixels_out<= 4*1024*1024;
   //bool force_st=false;
   if (nRunningThreads > 1 && !force_st) {
      for (int i = 0; i < nRunningThreads; i++) {
         params[i].setParams(n_pixels_out * i / nRunningThreads, n_pixels_out * (i + 1) / nRunningThreads);
         pool.addTask(threadEntry<T> , &params[i]);
      }

#ifdef DEBUG_TIME
      printf("nRunningThreads=%d; force_st=%d; #pxl_out=%ld\n",
             nRunningThreads, int(force_st), n_pixels_out);
      printf("before wait4all @  %.4f\n", read_timer());
#endif
      pool.wait4all();

#ifdef DEBUG_TIME
      printf("after wait4all @  %.4f\n", read_timer());
#endif
   }
   else {
      params[0].setParams(0, n_pixels_out);
#ifdef DEBUG_TIME
      printf("nRunningThreads=%d; force_st=%d; #pxl_out=%ld\n",
             nRunningThreads, int(force_st), n_pixels_out);
      printf("before call threadEntry directly @ %.4f\n", read_timer());
#endif
      threadEntry<T> (&params[0]);
#ifdef DEBUG_TIME
      printf("after call threadEntry directly @ %.4f\n", read_timer());
#endif
   }
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

  start_timer();

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables g and im\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  printf("before load var() @  %.4f\n", read_timer());
  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);
  printf("after load var() @  %.4f\n", read_timer());

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  printf("before save var() @  %.4f\n", read_timer());
  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);
  printf("after save var() @  %.4f\n", read_timer());

  return EXIT_SUCCESS;
}
#endif
