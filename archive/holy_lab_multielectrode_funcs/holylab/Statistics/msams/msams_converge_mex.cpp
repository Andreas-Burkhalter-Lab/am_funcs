#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#else
#include "mex.h"
#endif

//#include <math.h>
#include "msams_converge_core.cpp"

/* msams_converge_mex: move point(s) by adaptive meanshift, using landmarking
 *
 * Syntax:
 *   outdata = msams_converge_mex(x,lminfo,y,options)
 * where (syntax is similar to msams_converge)
 *   x is a d-by-N matrix of data points in d-dimensional space;
 *   lminfo is a landmark structure of the type returned by
 *     choose_landmarks;
 *   y is a d-by-q matrix of probe points in d-dimensional space;
 *   options is a structure which may have the following fields:
 *     min_to_check (default 3)
 *     factor (default 3)
 *     terminate_mode (default 'c')
 *     backtrack (default true)
 *     any_coordinate (default false)
 *     convergence_thresh (default 1e-12)
 *     max_iter (default 1000)
 *     n_threads (default # of cores)
 * and outdata is a structure which has the following fields:
 *   yf is d-by-q matrix containing the final position of the
 *     probe points
 *   closestDataIndex is the 1-by-q vector holding the index of the data
 *     point closest to the final position yf
 *   n is the 1-by-q vector containing the number of points
 *     contributing to each probe point at its final position
 *   R2 is the squared radius for each final probe point
 *   n_iter is the number of iterations required for convergence for
 *     each probe point
 *   convergedFlag is a flag indicating whether each point converged
 *     (0 = no convergence, 1 = converged, 2 = "converged" by cycling)
 *   settings is a copy of the options structure with defaults filled
 *     in (checking this is the most reliable way to learn what the
 *     defaults really are!)
 *   ntraj (present only for q = 1) is the sequence of # of contributing
 *     neighbors on each mean shift cycle
 *   ytraj (present only for q = 1) is a d-by-niterations matrix that
 *     stores the history of positions visited
 */

// Utility functions
//int is_scalar(const mxArray *m);
//double mcm_get_scalar_field(const mxArray *mxPtr,const char *name);
template <class T> 
void fillOptionalScalarField(const mxArray *m,const char *name,T *v);
template <class T>
void setFieldPtr(const mxArray *m,const char *name,T **d);
void mat2C(const mxArray *mxOptions,optionStruct &ops);
void C2mat(const optionStruct &ops,mxArray *mxOptions);
template <class Tdata,class Tint>
void mat2C(const mxArray *mxOut,outputStruct<Tdata,Tint> &out);


// Information specifying structures for Matlab's consumption
const char* settingsfields[] = {
  "min_to_check",
  "factor",
  "terminate_mode",
  "backtrack",
  "any_coordinate",
  "convergence_thresh",
  "max_iter",
  "n_threads"
};
const int n_settingsfields = 8;

const char* outputfields[] = {
  "yf",
  "closestDataIndex",
  "n",
  "R2",
  "n_iter",
  "convergedFlag",
  "settings",
  "ntraj",
  "ytraj"
};
const int n_outputfields = 7;  // add 2 if q==1 to get ntraj & ytraj

const int matlab_offset = 1;

// This is the "inner" matlab wrapper. It is called by mexFunction
// (below), but this one does all the real work.  It's templated so
// that it can work with a variety of data types.
template <class dataType>
void msams_wrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[],bool (*validationFunction)(const mxArray*),mxClassID allocationTypeFlag)
{
  dataType *x,*y;
  double *lmIndex;
  int N,d,q,i;
  const mxArray *curarg;
  mxArray *mxPtr, *mxPtr2, *mxOutput, *settingsStruct;
  landmarkStruct<dataType,double> lminfo;
  optionStruct ops;
  outputStruct<dataType,double> out;

  // Parse the inputs
  // Get the data points
  curarg = prhs[0];
  x = (dataType*) mxGetData(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);

  // Collect the landmark info
  curarg = prhs[1];
  if (!mxIsStruct(curarg))
    mexErrMsgTxt("msams_converge_mex: lminfo must be a structure");
  //Landmark info: landmarks field
  mxPtr = mxGetField(curarg,0,"landmarks");
  if (mxPtr == NULL)
    mexErrMsgTxt("msams_converge_mex: error parsing lminfo's 'landmarks' field");
  lminfo.n_landmarks = mxGetN(mxPtr);
  if (mxGetM(mxPtr) != d)
    mexErrMsgTxt("msams_converge_mex: landmarks must have the same dimensionality as the data");
  if (!validationFunction(mxPtr))
    mexErrMsgTxt("msams_converge_mex: landmarks must be of the same data type as the data points");
  lminfo.landmarks = (dataType*) mxGetData(mxPtr);
  //Landmark info: landmarkAssignment field
  mxPtr = mxGetField(curarg,0,"landmarkAssignment");
  if (mxPtr == NULL)
    mexErrMsgTxt("msams_converge_mex: error parsing lminfo's 'landmarkAssignment' field");
  if (mxGetNumberOfElements(mxPtr) != N)
    mexErrMsgTxt("msams_converge_mex: wrong number of entries in 'landmarkAssignment' field");
  lminfo.landmarkAssignment = mxGetPr(mxPtr);
  //Landmark info: radius_of_neighborhood field
  mxPtr = mxGetField(curarg,0,"radius_of_neighborhood");
  if (mxPtr == NULL)
    mexErrMsgTxt("msams_converge_mex: error parsing lminfo's 'radius_of_neighborhood' field");
  if (mxGetNumberOfElements(mxPtr) != lminfo.n_landmarks)
    mexErrMsgTxt("msams_converge_mex: wrong number of entries in 'radius_of_neighborhood' field");
  lminfo.landmarkR = mxGetPr(mxPtr);
  //Landmark info: landmarkList field
  mxPtr = mxGetField(curarg,0,"landmarkList");
  if (mxPtr == NULL)
    mexErrMsgTxt("msams_converge_mex: error parsing lminfo's 'landmarkList' field");
  if (!mxIsCell(mxPtr))
    mexErrMsgTxt("msams_converge_mex: lminfo's 'landmarkList' field must be a cell array");
  if (mxGetN(mxPtr) != lminfo.n_landmarks)
    mexErrMsgTxt("msams_converge_mex: landmarkList must have the same number of elements as there are landmarks");
  lminfo.landmarkList = (double **) mxMalloc(lminfo.n_landmarks * sizeof(double*));
  lminfo.n_landmarkList = (int *) mxMalloc(lminfo.n_landmarks * sizeof(int));
  for (i = 0; i < lminfo.n_landmarks; i++) {
    mxPtr2 = mxGetCell(mxPtr,i);
    lminfo.landmarkList[i] = mxGetPr(mxPtr2);
    lminfo.n_landmarkList[i] = mxGetNumberOfElements(mxPtr2);
  }
  lminfo.d = d;
  lminfo.index_offset = matlab_offset;

  // Get the probe points
  curarg = prhs[2];
  if (mxIsEmpty(curarg)) {
    // Passing an empty matrix uses the data points as the probe points
    y = x;
    q = N;
  } else {
    // Read the probe points
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !validationFunction(curarg))
      mexErrMsgTxt("msams_converge_mex: y must be of the same data type as the data points");
    if (mxGetM(curarg) != d)
      mexErrMsgTxt("msams_converge_mex: the number of rows in x and y must be the same");
    q = mxGetN(curarg);
    y = (dataType*) mxGetData(curarg);
  }

  /*
  landmarked_neighbors<dataType,double> lm_nbrs;
  lm_nbrs.allocate(lminfo.n_landmarks,N);
  lm_nbrs.initialize(y,x,lminfo);
  return;
  */

  // Parse the options
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsStruct(curarg))
      mexErrMsgTxt("msams_converge_mex: options must be a structure");
    mat2C(curarg,ops);
  }
  ops.n_min = (N >= ops.n_min) ? ops.n_min : N;  // set to min(n_min,N)
  ops.n_threads = (q < ops.n_threads) ? q : ops.n_threads;
  ops.n_threads = (ops.n_cpus < ops.n_threads) ? ops.n_cpus : ops.n_threads;

  // Print info
  /*
  mexPrintf("d %d, N %d, q %d, nlhs %d, nrhs %d\n",d,N,q,nlhs,nrhs);
  mexPrintf("factor %g, convergence_thresh %g\n",ops.factor,ops.convergence_thresh);
  mexPrintf("leapfrog %d, n_min %d, backtrack %d, any_coordinate %d, max_iter %d\n",ops.leapfrog,ops.n_min,ops.backtrack,ops.any_coordinate,ops.max_iter);
  mexPrintf("n_landmarks_in_list: ");
  for (i = 0; i < lminfo.n_landmarks; i++)
    mexPrintf("%d ",lminfo.n_landmarkList[i]);
  mexPrintf("\n");
  */

  // Set up the output; the fields depend on whether q == 1
  plhs[0] = mxCreateStructMatrix(1,1,n_outputfields+2*(q==1),outputfields);
  mxOutput = plhs[0];
  mxSetField(mxOutput,0,"yf",mxCreateNumericMatrix(d,q,allocationTypeFlag,mxREAL));
  mxSetField(mxOutput,0,"closestDataIndex",mxCreateDoubleMatrix(1,q,mxREAL));
  mxSetField(mxOutput,0,"n",mxCreateDoubleMatrix(1,q,mxREAL));
  mxSetField(mxOutput,0,"R2",mxCreateNumericMatrix(1,q,allocationTypeFlag,mxREAL));
  mxSetField(mxOutput,0,"n_iter",mxCreateDoubleMatrix(1,q,mxREAL));
  mxSetField(mxOutput,0,"convergedFlag",mxCreateDoubleMatrix(1,q,mxREAL));
  mat2C(mxOutput,out);
  memcpy(out.y,y,d*q*sizeof(dataType));   // copy the probe point positions
  settingsStruct = mxCreateStructMatrix(1,1,n_settingsfields,settingsfields);
  C2mat(ops,settingsStruct);
  mxSetField(mxOutput,0,"settings",settingsStruct);
  if (q == 1) {
    // Do ntraj and ytraj as temporary variables
    out.ntraj = (double *) mxMalloc(ops.max_iter * sizeof(double));
    out.xnbrI = (double *) mxMalloc(ops.max_iter * sizeof(double));
    out.ytraj = (dataType*) mxMalloc((ops.max_iter+1)*d*sizeof(dataType));
  }
  
  if (!msams_core(x,d,N,lminfo,q,out,ops))
    mexErrMsgTxt("msams_converge_mex: error creating threads");

  if (q == 1) {
    // Convert ntraj & ytraj to output mxArrays
    int n_iter = (int) out.n_iter[0];
    mxArray *ntrajMx,*ytrajMx;
    ntrajMx = mxCreateDoubleMatrix(1,n_iter,mxREAL);
    mxSetField(mxOutput,0,"ntraj",ntrajMx);
    ytrajMx = mxCreateNumericMatrix(d,n_iter+1,allocationTypeFlag,mxREAL);
    mxSetField(mxOutput,0,"ytraj",ytrajMx);
    memcpy(mxGetPr(ntrajMx),out.ntraj,n_iter*sizeof(double));
    memcpy(mxGetPr(ytrajMx),out.ytraj,(n_iter+1)*d*sizeof(dataType));
    mxFree(out.ntraj);
    mxFree(out.xnbrI);
    mxFree(out.ytraj);
  }

  mxFree(lminfo.landmarkList);
  mxFree(lminfo.n_landmarkList);
  return;
}

// The "outer" matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("msams_converge_mex: requires three or four inputs");
  if (nlhs != 1)
    mexErrMsgTxt("msams_converge_mex: requires one output");

  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("msams_converge_mex: x must be a real matrix");

  if (mxIsDouble(curarg))
    msams_wrapper<double>(nlhs,plhs,nrhs,prhs,&mxIsDouble,mxDOUBLE_CLASS);
  else if (mxIsSingle(curarg))
    msams_wrapper<float>(nlhs,plhs,nrhs,prhs,&mxIsSingle,mxSINGLE_CLASS);
  else
    mexErrMsgTxt("msams_converge_mex: x must be a single- or double-precision");
    
  return;
}

template <class T>
void fillOptionalScalarField(const mxArray *mxPtr,const char *name,T *v)
{
  const mxArray *fieldPtr;
 
  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr != NULL) {
    if (mxGetNumberOfElements(fieldPtr) != 1)
      mexErrMsgIdAndTxt("msams_converge_mex:field_parsing_error","msams_converge_mex: expect field '%s' to be a scalar",name);
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
  fillOptionalScalarField(mxOptions,"terminate_mode",&(ops.terminate_mode));
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
  const int sz = 1;
  mxArray *tm = mxCreateCharArray(1,&sz);
  mxChar *charData = (mxChar*) mxGetData(tm);
  charData[0] = ops.terminate_mode;
  mxSetField(mxOptions,0,"terminate_mode",tm);
  mxSetField(mxOptions,0,"backtrack",mxCreateScalarDouble(double(ops.backtrack)));
  mxSetField(mxOptions,0,"any_coordinate",mxCreateScalarDouble(double(ops.any_coordinate)));
  mxSetField(mxOptions,0,"convergence_thresh",mxCreateScalarDouble(ops.convergence_thresh));
  mxSetField(mxOptions,0,"max_iter",mxCreateScalarDouble(double(ops.max_iter)));
  mxSetField(mxOptions,0,"n_threads",mxCreateScalarDouble(double(ops.n_threads)));
}    

template <class Tdata,class Tint>
void mat2C(const mxArray *mxOut,outputStruct<Tdata,Tint> &out) {
  setFieldPtr(mxOut,"yf",&(out.y));
  setFieldPtr(mxOut,"closestDataIndex",&(out.closestDataIndex));
  setFieldPtr(mxOut,"n",&(out.n));
  setFieldPtr(mxOut,"R2",&(out.R2));
  setFieldPtr(mxOut,"n_iter",&(out.n_iter));
  setFieldPtr(mxOut,"convergedFlag",&(out.convergedFlag));
};


// For debugging & profiling: build a stand-alone application
#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 4;
  const int n_outputs = 1;
  mxArray *input[n_inputs];
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

  printf("Output file just before calling mexfcn: %s\n",fileout);
  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  // Save the outputs
  printf("Output file just before save: %s\n",fileout);
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
