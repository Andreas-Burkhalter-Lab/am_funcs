#include <unistd.h>  // for # of CPUs
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#else
#include "mex.h"
#endif

#include "msams_val_gradient_core.cpp"

// Copyright 2008 by Timothy E. Holy

/* The plan:
   [step,var,weights,gradratio] = msams_val_gradient(x,x0,kvec,max_threads)
 (max_threads is optional)
 Weights comes before gradratio because it takes extra computation to calculate the gradient, whereas the weights have to be computed anyway.
*/

typedef double dataType;

/*
template <class Tdata>
extern int msams_val_gradient_core(const Tdata *x,const Tdata *x0,const Tdata *kvec,int d,int N,Tdata *step,Tdata *var,Tdata *grad,Tdata *w,int n_threads);
*/

// The Matlab wrapper
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  dataType *x,*x0,*kvec;
  int d,N;
  int max_threads,ncpus;
  dataType *step,*var,*w,*grad;
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("msams_val_gradient: requires three or four inputs");
  if (nlhs < 2 || nlhs > 4)
    mexErrMsgTxt("msams_val_gradient: requires two to four outputs");

  // Parse the inputs
  // Get the data points
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("msams_val_gradient: x must be a real double-precision matrix");
  x = mxGetPr(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);
  if (N < 2)
    mexErrMsgTxt("msams_val_gradient: there must be at least two data points");

  // Get the probe point
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("msams_val_gradient: x0 must be a real double-precision vector");
  if (mxGetM(curarg) != d)
    mexErrMsgTxt("msams_val_gradient: the number of rows in x and x0 must be the same");
  if (mxGetN(curarg) != 1)
    mexErrMsgTxt("msams_val_gradient: handles only one probe point at a time");
  x0 = mxGetPr(curarg);

  // Get the metric coefficients
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("msams_val_gradient: kvec must be a real double-precision vector");
  if (mxGetNumberOfElements(curarg) != d)
    mexErrMsgTxt("msams_val_gradient: kvec must be of the same dimensionality as x and x0");
  kvec = mxGetPr(curarg);

  // Get max_threads
  ncpus = sysconf(_SC_NPROCESSORS_ONLN);
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
      mexErrMsgTxt("msams_val_gradient: max_threads must be a real integer");
    max_threads = (int) mxGetScalar(curarg);
  } else
    max_threads = ncpus;   // Use as many as there are CPUs
  // Don't use any more threads than CPUs
  if (max_threads > ncpus)
    max_threads = ncpus;


    
  // Set up storage for the outputs
  plhs[0] = mxCreateDoubleMatrix(d,1,mxREAL);
  step = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(d,1,mxREAL);
  var = mxGetPr(plhs[1]);
  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
    w = mxGetPr(plhs[2]);
  } else
    w = NULL;
  if (nlhs > 3) {
    plhs[3] = mxCreateDoubleMatrix(d,1,mxREAL);
    grad = mxGetPr(plhs[3]);
  } else
    grad = NULL;


  // Do the work!
  if (!msams_val_gradient_core(x,x0,kvec,d,N,step,var,grad,w,max_threads))
    mexErrMsgTxt("msams_val_gradient: the core function could not create threads");
}


/*
  

  if (q == 1) {
    // Do ntraj and ytraj as temporary variables
    out.ntraj = (double *) mxMalloc(ops.max_iter * sizeof(double));
    out.ytraj = (double *) mxMalloc((ops.max_iter+1)*d*sizeof(double));
  }
  
  if (!msams_core(x,d,N,lminfo,q,out,ops))
    mexErrMsgTxt("msams_val_gradient: error creating threads");

  if (q == 1) {
    // Convert ntraj & ytraj to output mxArrays
    int n_iter = (int) out.n_iter[0];
    mxArray *ntrajMx,*ytrajMx;
    ntrajMx = mxCreateDoubleMatrix(1,n_iter,mxREAL);
    mxSetField(mxOutput,0,"ntraj",ntrajMx);
    ytrajMx = mxCreateDoubleMatrix(d,n_iter+1,mxREAL);
    mxSetField(mxOutput,0,"ytraj",ytrajMx);
    memcpy(mxGetPr(ntrajMx),out.ntraj,n_iter*sizeof(double));
    memcpy(mxGetPr(ytrajMx),out.ytraj,(n_iter+1)*d*sizeof(double));
    mxFree(out.ntraj);
    mxFree(out.ytraj);
  }

  mxFree(lminfo.landmarkList);
  mxFree(lminfo.n_landmarkList);
  return;
}
*/

/*
template <class T>
void fillOptionalScalarField(const mxArray *mxPtr,const char *name,T *v)
{
  const mxArray *fieldPtr;
 
  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr != NULL) {
    if (mxGetNumberOfElements(fieldPtr) != 1)
      mexErrMsgIdAndTxt("msams_val_gradient:field_parsing_error","msams_val_gradient: expect field '%s' to be a scalar",name);
    *v = (T) mxGetScalar(fieldPtr);
    //mexPrintf("Field %s was set to %g\n",name,mxGetScalar(fieldPtr));
  }
}
  
template <class T>
void setFieldPtr(const mxArray *mxPtr,const char *name,T **d)
{
  const mxArray *fieldPtr;
 
  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr != NULL)
    *d = (T*) mxGetData(fieldPtr);
}

void mat2C(const mxArray *mxOptions,optionStruct &ops)
{
  fillOptionalScalarField(mxOptions,"min_to_check",&(ops.n_min));
  fillOptionalScalarField(mxOptions,"factor",&(ops.factor));
  fillOptionalScalarField(mxOptions,"leapfrog",&(ops.leapfrog));
  fillOptionalScalarField(mxOptions,"backtrack",&(ops.backtrack));
  fillOptionalScalarField(mxOptions,"any_coordinate",&(ops.any_coordinate));
  fillOptionalScalarField(mxOptions,"convergence_thresh",&(ops.convergence_thresh));
  fillOptionalScalarField(mxOptions,"max_iter",&(ops.max_iter));
  fillOptionalScalarField(mxOptions,"n_threads",&(ops.n_threads));
}

void C2mat(const optionStruct &ops,mxArray *mxOptions)
{
  mxSetField(mxOptions,0,"min_to_check",mxCreateScalarDouble(double(ops.n_min)));
  mxSetField(mxOptions,0,"factor",mxCreateScalarDouble(ops.factor));
  mxSetField(mxOptions,0,"leapfrog",mxCreateScalarDouble(double(ops.leapfrog)));
  mxSetField(mxOptions,0,"backtrack",mxCreateScalarDouble(double(ops.backtrack)));
  mxSetField(mxOptions,0,"any_coordinate",mxCreateScalarDouble(double(ops.any_coordinate)));
  mxSetField(mxOptions,0,"convergence_thresh",mxCreateScalarDouble(ops.convergence_thresh));
  mxSetField(mxOptions,0,"max_iter",mxCreateScalarDouble(double(ops.max_iter)));
  mxSetField(mxOptions,0,"n_threads",mxCreateScalarDouble(double(ops.n_threads)));
}    

template <class Tdata,class Tint>
void mat2C(const mxArray *mxOut,outputStruct<Tdata,Tint> &out) {
  setFieldPtr(mxOut,"yf",&(out.y));
  setFieldPtr(mxOut,"closestLandmark",&(out.closestLandmark));
  setFieldPtr(mxOut,"n",&(out.n));
  setFieldPtr(mxOut,"R2",&(out.R2));
  setFieldPtr(mxOut,"n_iter",&(out.n_iter));
  setFieldPtr(mxOut,"convergedFlag",&(out.convergedFlag));
};
*/

// For debugging & profiling: build a stand-alone application
#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 3;
  const int n_outputs = 1;
  const mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "x",
    "lminfo",
    "y",
    "ops"
  };
  const char *output_names[] = {
    "out"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables x, y, lminfo, and ops\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  //printf("Output file just before calling mexfcn: %s\n",fileout);
  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,input);

  // Save the outputs
  //printf("Output file just before save: %s\n",fileout);
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
