#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <iostream>
#include <valgrind/callgrind.h>
#else
#include "mex.h"
#endif

#include "Eigen/Geometry"

using namespace Eigen;

//
// Syntax: see rotate_quaternion.m
//

/*
 * This function is commented out because it is more efficient (when
 * rotating multiple vectors) to first convert the quaternion to a
 * rotation matrix. (The main advantage of quaternions comes in
 * composing rotations.)
 *
template <class dataType>
void rotate_quaternion(const dataType* X,int N,const Quaternion<dataType> &q,dataType *Xr)
{
  Quaternion<dataType> qc = q.conjugate();
  Quaternion<dataType> x,xr;
  const dataType *Xend = X + 3*N;

  for (; X < Xend; X += 3, Xr += 3) {
    x.w() = 0;
    x.x() = X[0];
    x.y() = X[1];
    x.z() = X[2];
    xr = q*x*qc;
    Xr[0] = xr.x();
    Xr[1] = xr.y();
    Xr[2] = xr.z();
  }
}
*/
  
// This is the "inner" matlab wrapper. It is called by mexFunction
// (below), but this one does all the real work.  It's templated so
// that it can work with a variety of data types.
template <class dataType>
void mexWrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[],mxClassID dataClassID)
{
  const mxArray *curarg;

  // Define a collection of 3-vectors using Matlab memory ordering
  typedef Matrix<dataType, 3, Dynamic, ColMajor> V3;
  typedef const Matrix<dataType, 3, Dynamic, ColMajor> V3const;

#ifdef MAIN
  CALLGRIND_START_INSTRUMENTATION;
#endif

  // Parse the inputs
  // X
  curarg = prhs[0];
  if (mxGetM(curarg) != 3)
    mexErrMsgTxt("X must have 3 rows");
  const dataType *Xp = (dataType*) mxGetData(curarg);
  // Convert to Eigen representation
  int N = mxGetN(curarg);
  Map<V3const> X(Xp,3,N);

  // q
  curarg = prhs[1];
  if (mxGetClassID(curarg) != dataClassID)
    mexErrMsgTxt("q must have the same data type as X");
  if (mxGetNumberOfElements(curarg) != 4)
    mexErrMsgTxt("q must have 4 elements");
  const dataType *qp = (dataType*) mxGetData(curarg);
  Quaternion<dataType> q(qp[0],qp[1],qp[2],qp[3]);
  // It's more efficient to first convert to a rotation matrix than it
  // is to use quaternion multiplication on each input vector.
  Matrix<dataType, 3, 3> R = q.toRotationMatrix();

  // Allocate storage for the output
  plhs[0] = mxCreateNumericMatrix(3,N,dataClassID,mxREAL);
  dataType *Xrp = (dataType *) mxGetData(plhs[0]);
  Map<V3> Xr(Xrp,3,N);

  // Perform the calculation
  Xr = R*X;

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
    mexErrMsgTxt("Requires 2 inputs, X and q");

  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("X must be a real numeric array");

  if (mxIsDouble(curarg))
    mexWrapper<double>(nlhs,plhs,nrhs,prhs,mxDOUBLE_CLASS);
  else if (mxIsSingle(curarg))
    mexWrapper<float>(nlhs,plhs,nrhs,prhs,mxSINGLE_CLASS);
  else
    mexErrMsgTxt("X must be a single- or double-precision");
    
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
