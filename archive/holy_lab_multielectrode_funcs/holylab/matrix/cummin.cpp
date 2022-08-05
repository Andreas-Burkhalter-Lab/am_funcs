#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <valgrind/callgrind.h>
#else
#include "mex.h"
#endif

#include "imiterators.cxx"

// Compile with mex -I../images/ cummin.cpp

//
// Syntax: cm = cummin(X,dim)
//

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
  // X
  curarg = prhs[0];
  const dataType *Xp = (dataType*) mxGetData(curarg);
  int n_dims = mxGetNumberOfDimensions(curarg);
  const mwSize *sz = mxGetDimensions(curarg);

  // dim
  int dim = 1;
  if (nrhs > 1) {
    curarg = prhs[1];
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("dim must be a scalar");
    dim = (int) mxGetScalar(curarg);
  }

  // Allocate storage for the output
  plhs[0] = mxCreateNumericArray(n_dims,sz,dataClassID,mxREAL);
  dataType* cmp = (dataType*) mxGetData(plhs[0]);

  // Perform the calculation
  pixIterator pI(sz,n_dims);
  dim--;
  for ( ; !pI.at_end(); pI++)
    if (pI.coord(dim) == 0)
      cmp[pI] = Xp[pI];
    else {
      dataType tmp1 = Xp[pI];
      dataType tmp2 = cmp[pI-pI.dimSkip(dim)];
      cmp[pI] = (tmp1 < tmp2) ? tmp1 : tmp2;
    }
    
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
  if (nrhs > 2 || nrhs < 1)
    mexErrMsgTxt("Requires 1-2 inputs, X and dim");

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
    "dim"
  };
  const char *output_names[] = {
    "cm"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables X and dim\n",argv[0]);
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
