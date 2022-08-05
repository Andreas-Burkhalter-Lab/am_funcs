// Define MAIN while compiling to create a stand-alone command-line
// program.
// To compile as a MEX file:
//    mex -I.. imflow_dotproduct_mex.cpp ../image_utilities.cpp
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#endif

#include "image_utilities.h"

#define UNIT_OFFSET 1

/* Syntax:
 *   map = imflow_dotproduct(im)
 */

// forward declarations
void imflow_work_double_int(const double *im, int n_values, int n_dims,const int *imsz,int *map);
void imflow_work_double_double(const double *im, int n_values, int n_dims,const int *imsz,double *map);
void imflow_work_float_int(const float *im, int n_values, int n_dims,const int *imsz,int *map);
void imflow_work_float_double(const float *im, int n_values, int n_dims,const int *imsz,double *map);


/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int n_values,n_dims,n_dims_orig;
  const int *imsz_orig;  // the size of im
  int *imsz;     // the size of im, without any unity dimensions
  const mxArray *curarg;
  bool isdouble,map_is_int;
  long n_pixels;
  void *im, *map;
  int i;

  if (nrhs != 1)
    mexErrMsgTxt("imflow_dotproduct_mex: requires one input, and the first dimension is interpreted as a value dimension");
  if (nlhs != 1)
    mexErrMsgTxt("imflow_dotproduct_mex: requires one output");

  // Parse the input
  // im
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("imflow_dotproduct_mex: A must be a real array");
  if (mxIsSingle(curarg))
    isdouble = false;
  else if (mxIsDouble(curarg))
    isdouble = true;
  else
    mexErrMsgTxt("imflow_dotproduct_mex: A must be single or double");
  if (mxIsEmpty(curarg))
    mexErrMsgTxt("imflow_dotproduct_mex: does not work on empty arrays");

  im = mxGetData(curarg);
  n_dims_orig = mxGetNumberOfDimensions(curarg);
  imsz_orig = mxGetDimensions(curarg);
  // First dimension must be value dimension
  n_values = imsz_orig[0];
  // Now for the spatial part just skip over the first dimension
  imsz_orig++;
  n_dims_orig--; 
  imsz = (int *) mxMalloc(n_dims_orig*sizeof(int));
  n_dims = skip_unity_dimensions(imsz_orig,n_dims_orig,imsz,n_dims_orig);
  n_pixels = 1;
  for (i = 0; i < n_dims; i++)
    n_pixels *= imsz[i];

  /*
  // n_threads
  n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  //n_threads = n_cpus;
  n_threads = 1;  // Currently not much benefit to multiple threads
  if (nrhs > 2) {
    curarg = prhs[2];
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("imflow_dotproduct_mex: n_threads must be a scalar");
    n_threads = (int) mxGetScalar(curarg);
    if (n_threads < 1)
      n_threads = 1;
    if (n_threads > n_cpus)
      n_threads = n_cpus;
  }
  */

  // Set up the output and do the work
  if (n_pixels < 2e9) {
    map_is_int = true;
    plhs[0] = mxCreateNumericArray(n_dims_orig,imsz_orig,
				 mxINT32_CLASS,mxREAL);
  } else {
    map_is_int = false;
    plhs[0] = mxCreateNumericArray(n_dims_orig,imsz_orig,
				 mxDOUBLE_CLASS,mxREAL);
  }
  map = mxGetData(plhs[0]);


  if (isdouble && map_is_int)
    imflow_work_double_int((double *) im,n_values,n_dims,imsz,(int *) map);
  else if (isdouble)
    imflow_work_double_double((double *) im,n_values,n_dims,imsz,(double *) map);
  else if (map_is_int)
    imflow_work_float_int((float *) im,n_values,n_dims,imsz,(int *) map);
  else
    imflow_work_float_double((float *) im,n_values,n_dims,imsz,(double *) map);

  return;
}


template <class Tdata,class Tmap>
void imflow_work(const Tdata *im,int n_values,int n_dims,const int *imsz,Tmap *map)
{
  long i,n_pixels;
  int dimIndex;
  int j;
  int n_nbrs;
  int *coords;
  long *dimSkip;
  long *offsets;
  long maxvalOffset;
  Tdata maxval,dp;
  const Tdata *pPix,*thisPixel,*thisPixelEnd,*nbrPixel;
  long imin,imax;

  //
  // Initialize indexing and offsets
  //
  // First calculate the # of nearest-neighbors and the memory offset
  // along each dimension
  // Neighbors are all points inside a 3-by-3-by... rectangular block
  // centered on the current pixel
  n_nbrs = 1;
  coords = new int[n_dims+1];  // +1 as a sentinel for increment-with-carry
  dimSkip = new long[n_dims];  // stores the memory increment along each axis
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    n_nbrs *= 3;
    coords[dimIndex] = -1;  // start all coords at -1
    if (dimIndex == 0)
      dimSkip[dimIndex] = 1;
    else
      dimSkip[dimIndex] = dimSkip[dimIndex-1] * imsz[dimIndex-1];
  }
  n_pixels = dimSkip[n_dims-1] * imsz[n_dims-1];

  // Calculate the pixel index offset corresponding to each neighbor
  // Note that memory offsets will be n_values * offsets[i]
  offsets = new long[n_nbrs];
  for (i = 0; i < n_nbrs; i++) {
    // Calculate the pixel index offset
    offsets[i] = 0;
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      offsets[i] += coords[dimIndex] * dimSkip[dimIndex];
    // Increment the coords to the next neighbor
    coords[0]++;
    dimIndex = 0;
    while (coords[dimIndex] > 1 && dimIndex < n_dims) {
      coords[dimIndex] = -1;
      dimIndex++;
      coords[dimIndex]++;
    }
  }
    
  // Iterate over all pixels of the image. It will be faster if we
  // don't worry about image boundaries (fix those later with
  // killedges)
  imin = -offsets[0];  //  offsets[0] is all coords = -1
  imax = n_pixels - offsets[n_nbrs-1]; // offsets[n_nbrs-1] is all coords = 1
  for (i = imin; i < imax; i++) {  // iteration over "memory-valid" pixels
    thisPixel = im + i*n_values;
    thisPixelEnd = thisPixel + n_values;
    // Calculate dot product with self Since this comes first and
    // others must exceed this value, ties will result in self-mapping
    maxvalOffset = 0;  // by default map to self
    dp = 0;
    for (pPix = thisPixel; pPix < thisPixelEnd; pPix++)
      dp += *pPix * *pPix;
    maxval = dp;
    // Calculate dot product with others
    for (j = 0; j < n_nbrs; j++) {
      nbrPixel = im + (i+offsets[j])*n_values;
      dp = 0;
      for (pPix = thisPixel; pPix < thisPixelEnd; pPix++,nbrPixel++)
	dp += *pPix * *nbrPixel;
      if (dp > maxval) {
	maxvalOffset = offsets[j];
	maxval = dp;
      }
    }
    map[i] = Tmap(i+UNIT_OFFSET+maxvalOffset);
  }

  delete[] offsets;
  delete[] dimSkip;
  delete[] coords;
}

void imflow_work_double_int(const double *im,int n_values,int n_dims,const int *imsz,int *map)
{
  imflow_work<double,int>(im,n_values,n_dims,imsz,map);
}

void imflow_work_double_double(const double *im,int n_values,int n_dims,const int *imsz,double *map)
{
  imflow_work<double,double>(im,n_values,n_dims,imsz,map);
}

void imflow_work_float_int(const float *im,int n_values,int n_dims,const int *imsz,int *map)
{
  imflow_work<float,int>(im,n_values,n_dims,imsz,map);
}

void imflow_work_float_double(const float *im,int n_values,int n_dims,const int *imsz,double *map)
{
  imflow_work<float,double>(im,n_values,n_dims,imsz,map);
}




#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 3;
  const int n_outputs = 1;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "A",
    "dimFlag",
    "n_threads"
  };
  const char *output_names[] = {
    "Aout"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables A, dimFlag, and n_threads\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  printf("Output file just before calling mexfcn: %s\n",fileout);
  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  // Save the outputs
  printf("Output file just before save: %s\n",fileout);
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
