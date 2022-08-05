#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#else
#include "mex.h"
#endif

#include "Eigen/Geometry"

using namespace Eigen;

//
// Syntax: see multiply_quaternion.m
//

template <class dataType>
void multiply_quaternion(const dataType* a,const dataType *b,int N,dataType *ab)
{
  Quaternion<dataType> qa, qb, qab;
  const dataType *aEnd = a + 4*N;

  for (; a < aEnd; a += 4, b += 4, ab += 4) {
    qa.w() = a[0];
    qa.x() = a[1];
    qa.y() = a[2];
    qa.z() = a[3];
    qb.w() = b[0];
    qb.x() = b[1];
    qb.y() = b[2];
    qb.z() = b[3];
    qab = qa*qb;
    ab[0] = qab.w();
    ab[1] = qab.x();
    ab[2] = qab.y();
    ab[3] = qab.z();
  }
}
  
// This is the "inner" matlab wrapper. It is called by mexFunction
// (below), but this one does all the real work.  It's templated so
// that it can work with a variety of data types.
template <class dataType>
void mexWrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[],mxClassID dataClassID)
{
  const mxArray *curarg;

#ifdef MAIN
  CALLGRIND_START_INSTRUMENTATION;
#endif

  // Parse the inputs
  // a
  curarg = prhs[0];
  if (mxGetM(curarg) != 4)
    mexErrMsgTxt("a must have 4 rows");
  const dataType *a = (dataType*) mxGetData(curarg);
  int N = mxGetN(curarg);

  // b
  curarg = prhs[1];
  if (mxGetClassID(curarg) != dataClassID)
    mexErrMsgTxt("b must have the same data type as a");
  if (mxGetM(curarg) != 4)
    mexErrMsgTxt("b must have 4 rows");
  if (mxGetN(curarg) != N)
    mexErrMsgTxt("b must have the same number of quaternions as a");
  const dataType *b = (dataType*) mxGetData(curarg);

  // Allocate storage for the output
  plhs[0] = mxCreateNumericMatrix(4,N,dataClassID,mxREAL);
  dataType *ab = (dataType *) mxGetData(plhs[0]);

  // Perform the calculation
  multiply_quaternion<dataType>(a,b,N,ab);

#ifdef MAIN
  CALLGRIND_STOP_INSTRUMENTATION;
#endif
}

// The "outer" matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nlhs > 1)
    mexErrMsgTxt("Only one output argument is provided");
  if (nrhs != 2)
    mexErrMsgTxt("Requires 2 inputs, a and b");

  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("a must be a real numeric array");

  if (mxIsDouble(curarg))
    mexWrapper<double>(nlhs,plhs,nrhs,prhs,mxDOUBLE_CLASS);
  else if (mxIsSingle(curarg))
    mexWrapper<float>(nlhs,plhs,nrhs,prhs,mxSINGLE_CLASS);
  else
    mexErrMsgTxt("a must be a single- or double-precision");
    
  return;
}


// For debugging & profiling: build a stand-alone application
#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 2;
  const int n_outputs = 1;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "X",
    "q"
  };
  const char *output_names[] = {
    "Xr"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables dX, p, options, and sortOrder\n",argv[0]);
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
