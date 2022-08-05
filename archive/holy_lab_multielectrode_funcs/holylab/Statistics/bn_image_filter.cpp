#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <stdlib.h>
#else
#include "mex.h"
#endif

#include "imiterators.cxx"
//#include "newdelete.h"
#include <math.h>    // for round()
#include <unistd.h>  // for # of CPUs
#include <vector>

using namespace std;

// bn_image_filter: filter an image using balanced neighborhood criterion
//
// Syntax: see bn_image_filter.m
//

struct optionStruct {
  double z;
  int n_threads;
  int n_cpus;

  optionStruct() {
    z = 3;
    n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    n_threads = n_cpus;
    if (n_threads > 10)
      n_threads = 10;
  }
};

// Utility functions
//int is_scalar(const mxArray *m);
//double mcm_get_scalar_field(const mxArray *mxPtr,const char *name);
template <class T> 
void fillOptionalScalarField(const mxArray *m,const char *name,T *v);
void mat2C(const mxArray *mxOptions,optionStruct &ops);
void C2mat(const optionStruct &ops,mxArray *mxOptions);
template <class dataType>
void bn_image_filter_work(int n_dims,int n_spatial_dims,const int *szint,const dataType *imdata,const dataType *imfilteredIn,const vector<long> &memoffset,const vector<bool> &groupEnd,const vector<int> &n_lookup,const vector<int> &left,const vector<int> &right,const optionStruct &ops,double *n,dataType *imfilteredOut,double *n_triggered,vector<dataType> &dVsum);
void initialize_nlookup(const vector<bool> &groupEnd,const optionStruct &ops,vector<int> &n_lookup);
void initialize_edgetest(const vector<int> &coordoffset,const vector<bool> &groupEnd,const int *sz,int n_spatial_dims,vector<int> &left,vector<int> &right);



// Information specifying structures for interfacing with Matlab
const char* settingsfields[] = {
  "z",
  "n_threads"
};
const int n_settingsfields = 2;

//const int matlab_offset = 1;

// This is the "inner" matlab wrapper. It is called by mexFunction
// (below), but this one does all the real work.  It's templated so
// that it can work with a variety of data types.
template <class dataType>
void msams_wrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[],mxClassID allocationTypeFlag)
{
  const mxArray *curarg,*memoffsetMx,*coordoffsetMx;
  const dataType *imdata,*imfilteredIn;
  dataType *imfilteredOut;
  double *n,*n_triggered;
  const mwSize *sz,*sztmp;
  int n_dims,n_spatial_dims,n_value_dims,n_values;
  optionStruct ops;
  int i,j,dimIndex,curarg_n,n_groups,n_in_group;

  // Parse the inputs
  // Get the base image
  curarg = prhs[0];
  imdata = (dataType*) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz = mxGetDimensions(curarg);

  // Get the filtered image, if supplied
  curarg_n = 1;
  curarg = prhs[curarg_n];
  if (!mxIsNumeric(curarg))
    imfilteredIn = imdata;
  else {
    imfilteredIn = (dataType *) mxGetData(curarg);
    curarg_n++;
    if (mxGetNumberOfDimensions(curarg) != n_dims)
      mexErrMsgTxt("# of dimensions in filtered image does not match imdata");
    sztmp = mxGetDimensions(curarg);
    for (i = 0; i < n_dims; i++)
      if (sz[i] != sztmp[i])
	mexErrMsgTxt("Dimensionality of imfilteredIn does not match imdata");
  }

  // memoffset and coordoffset
  vector<long> memoffset;
  vector<bool> groupEnd;
  vector<int> coordoffset;
  memoffsetMx = prhs[curarg_n];
  coordoffsetMx = prhs[curarg_n+1];
  curarg_n += 2;
  if (!mxIsCell(memoffsetMx) || !mxIsCell(coordoffsetMx))
    mexErrMsgTxt("memoffset and coordoffset should be a cell arrays");
  n_groups = mxGetNumberOfElements(memoffsetMx);
  if (mxGetNumberOfElements(coordoffsetMx) != n_groups)
    mexErrMsgTxt("memoffset and coordoffset must have the same size");
  for (i = 0; i < n_groups; i++) {
    const mxArray *thisElementMx = mxGetCell(memoffsetMx,i);
    if (!mxIsDouble(thisElementMx))
      mexErrMsgTxt("Entries in memoffset must be of type double");
    n_in_group = mxGetNumberOfElements(thisElementMx);
    double *thisElement = mxGetPr(thisElementMx);
    for (j = 0; j < n_in_group; j++) {
      memoffset.push_back(thisElement[j]);
      groupEnd.push_back(j == n_in_group-1);
    }
    thisElementMx = mxGetCell(coordoffsetMx,i);
    if (!mxIsDouble(thisElementMx))
      mexErrMsgTxt("Entries in coordoffset must be of type double");
    if (mxGetN(thisElementMx) != n_in_group)
      mexErrMsgTxt("Entries in coordoffset must have the same number of columns as in memoffset");
    if (i == 0)
      n_spatial_dims = mxGetM(thisElementMx);
    else
      if (mxGetM(thisElementMx) != n_spatial_dims)
	mexErrMsgTxt("Entries in coordoffset must have a consistent number of dimensions");
    thisElement = mxGetPr(thisElementMx);
    for (j = 0; j < n_in_group; j++)
      for (dimIndex = 0; dimIndex < n_spatial_dims; dimIndex++,thisElement++)
	coordoffset.push_back(*thisElement);
  }
  n_value_dims = n_dims-n_spatial_dims;
  n_values = 1;
  for (dimIndex = 0; dimIndex < n_value_dims; dimIndex++)
    n_values *= sz[dimIndex];
  
  // options
  if (curarg_n < nrhs) {
    curarg = prhs[curarg_n];
    if (!mxIsStruct(curarg))
      mexErrMsgTxt("options must be a structure");
    mat2C(curarg,ops);
  }
  if (ops.n_threads > ops.n_cpus)
    ops.n_threads = ops.n_cpus;
  
  // Allocate outputs
  if (nlhs > 0) {
    // n
    plhs[0] = mxCreateNumericArray(n_spatial_dims,sz+n_value_dims,mxDOUBLE_CLASS,mxREAL);
    n = mxGetPr(plhs[0]);
  } else
    return;
  if (nlhs > 1) {
    // imfilteredOut
    plhs[1] = mxCreateNumericArray(n_dims,sz,allocationTypeFlag,mxREAL);
    imfilteredOut = (dataType *) mxGetData(plhs[1]);
  } else
    imfilteredOut = NULL;
  if (nlhs > 2) {
    // n_triggered
    plhs[2] = mxCreateNumericArray(n_spatial_dims,sz+n_value_dims,mxDOUBLE_CLASS,mxREAL);
    n_triggered = mxGetPr(plhs[2]);
  } else
    n_triggered = NULL;
  // Allocate temporary storage
  vector<dataType> dVsum(n_values);

  // Put sz info into an integer type
  int *szint = (int*) mxMalloc(sizeof(int)*n_dims);
  for (i = 0; i < n_dims; i++)
    szint[i] = sz[i];

  // Initialize n_lookup, the lookup function for the backtracking criterion
  vector<int> n_lookup;
  initialize_nlookup(groupEnd,ops,n_lookup);
  // Initialize left & right, the tests for coordinate values going
  // over-the-edge
  vector<int> left,right;
  initialize_edgetest(coordoffset,groupEnd,szint+n_value_dims,n_spatial_dims,left,right);

  // Print info
  /*
  mexPrintf("n_dims %d, n_spatial_dims %d, n_lookup: ",n_dims,n_spatial_dims);
  vector<int>::iterator i1;
  for (i1 = n_lookup.begin(); i1 != n_lookup.end(); i1++)
    mexPrintf("%d ",*i1);
  mexPrintf("\n");
  vector<int>::iterator ii;
  mexPrintf("left:\n");
  for (dimIndex = 1, ii = left.begin(); ii != left.end(); ii++,dimIndex++) {
    mexPrintf("%d ",*ii);
    if (dimIndex == n_spatial_dims) {
      dimIndex = 0;
      mexPrintf("\n");
    }
  }
  mexPrintf("right:\n");
  for (dimIndex = 1, ii = right.begin(); ii != right.end(); ii++,dimIndex++) {
    mexPrintf("%d ",*ii);
    if (dimIndex == n_spatial_dims) {
      dimIndex = 0;
      mexPrintf("\n");
    }
  }
  */

  // Call the work function
  bn_image_filter_work(n_dims,n_spatial_dims,szint,imdata,imfilteredIn,memoffset,groupEnd,n_lookup,left,right,ops,n,imfilteredOut,n_triggered,dVsum);

  // Free temporaries
  mxFree(szint);
}

// The "outer" matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 5)
    mexErrMsgTxt("Requires three to five inputs");

  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("imdata must be a real numeric array");

  if (mxIsDouble(curarg))
    msams_wrapper<double>(nlhs,plhs,nrhs,prhs,mxDOUBLE_CLASS);
  else if (mxIsSingle(curarg))
    msams_wrapper<float>(nlhs,plhs,nrhs,prhs,mxSINGLE_CLASS);
  else
    mexErrMsgTxt("imdata must be a single- or double-precision");
    
  return;
}

template <class T>
void fillOptionalScalarField(const mxArray *mxPtr,const char *name,T *v)
{
  const mxArray *fieldPtr;
 
  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr != NULL) {
    if (mxGetNumberOfElements(fieldPtr) != 1)
      mexErrMsgIdAndTxt("bn_image_filter:field_parsing_error","bn_image_filter: expect field '%s' to be a scalar",name);
    *v = (T) mxGetScalar(fieldPtr);
    //mexPrintf("Field %s was set to %g\n",name,mxGetScalar(fieldPtr));
  }
}
  
void mat2C(const mxArray *mxOptions,optionStruct &ops)
{
  fillOptionalScalarField(mxOptions,"z",&(ops.z));
  fillOptionalScalarField(mxOptions,"n_threads",&(ops.n_threads));
}

void C2mat(const optionStruct &ops,mxArray *mxOptions)
{
  mxSetField(mxOptions,0,"z",mxCreateScalarDouble(double(ops.z)));
  mxSetField(mxOptions,0,"n_threads",mxCreateScalarDouble(double(ops.n_threads)));
}    

void initialize_nlookup(const vector<bool> &groupEnd,const optionStruct &ops,vector<int> &n_lookup)
{
  vector<int> n_cum;  // will hold the cumulative # of nbrs up to given group#
  int i,n_target,dist,mindist,nbest;
  double z2;
  
  z2 = ops.z*ops.z;

  n_lookup.clear();
  for (i = 0; i < groupEnd.size(); i++)
    if (groupEnd[i])
      n_cum.push_back(i+1);
  vector<int>::iterator i1,i2;
  for (i1 = n_cum.begin(); i1 != n_cum.end(); i1++) {
    n_target = (int) round((*i1-z2)/2);
    // Find the value of n_cum closest to n_target
    mindist = n_cum[n_cum.size()-1]+(int) z2;// sentinel value, bigger than all
    for (i2 = n_cum.begin(); i2 != n_cum.end(); i2++) {
      dist = *i2-n_target;
      if (dist < 0)
	dist = -dist;  // abs(dist)
      if (dist < mindist) {
	mindist = dist;
	nbest = *i2;
      }
    }
    n_lookup.push_back(nbest);
  }
}
	  
void initialize_edgetest(const vector<int> &coordoffset,const vector<bool> &groupEnd,const int *sz,int n_spatial_dims,vector<int> &left,vector<int> &right)
{
  vector<int>::const_iterator ci;
  vector<bool>::const_iterator gi;
  vector<int>::iterator li,litmp;
  vector<int>::iterator ri,ritmp;
  bool startFlag;
  int dimIndex,n_groups;

  for (gi = groupEnd.begin(),n_groups=0; gi != groupEnd.end(); gi++)
    if (*gi)
      n_groups++;
  left.clear();
  right.clear();
  left.reserve(n_groups*n_spatial_dims);
  right.reserve(n_groups*n_spatial_dims);

  ci = coordoffset.begin();
  li = left.begin();
  ri = right.begin();
  startFlag = true;

  for (gi = groupEnd.begin(); gi != groupEnd.end(); gi++) {
    if (startFlag) {
      // Put sentinel values into left and right
      for (dimIndex = 0; dimIndex < n_spatial_dims; dimIndex++) {
	left.push_back(0);
	right.push_back(sz[dimIndex]);
      }
      startFlag = false;
      if (&li[0] == NULL) {
	li = left.begin();
	ri = right.begin();
      }
    }
    // Determine whether the current coordoffset exceeds the
    // boundaries established for this group so far
    for (dimIndex=0,litmp=li,ritmp=ri; dimIndex < n_spatial_dims; dimIndex++,ci++,litmp++,ritmp++) {
      int thisco = *ci;
      if (thisco + *litmp < 0)
	*litmp=-thisco;
      if (thisco + *ritmp > sz[dimIndex])
	*ritmp=sz[dimIndex]-thisco;
    }
    // If we are at the end of a group, prepare for the next group
    if (*gi) {
      li += n_spatial_dims;
      ri += n_spatial_dims;
      startFlag = true;
    }
  }
}    
      
template <typename dataType>
void bn_image_filter_work(int n_dims,int n_spatial_dims,const int *sz,const dataType *imdata,const dataType *imfilteredIn,const vector<long> &memoffset,const vector<bool> &groupEnd,const vector<int> &n_lookup,const vector<int> &left,const vector<int> &right,const optionStruct &ops,double *n,dataType *imfilteredOut,double *n_triggered,vector<dataType> &dVsum)
{
  int dimIndex;
  const dataType *imdatatmp;
  dataType dVsum2;
  typedef typename vector<dataType>::iterator dataTypeIterator;

  // Initialization
  int n_value_dims = n_dims-n_spatial_dims;
  long n_values = 1;
  for (dimIndex = 0; dimIndex < n_value_dims; dimIndex++)
    n_values *= sz[dimIndex];
  double z2 = ops.z*ops.z;

  pixIterator pI(sz+n_value_dims,n_spatial_dims);

  for (; !pI.at_end(); pI++,imdata+=n_values,imfilteredIn+=n_values,n++) {
    // Initialize for current pixel
    //mexPrintf("\nPixel ");
    //for (dimIndex = 0; dimIndex < n_spatial_dims; dimIndex++)
    //  mexPrintf("%d ",pI.coord(dimIndex));
    //mexPrintf("\n(Memory offsets) values: ");
    *n = 0;
    fill(dVsum.begin(),dVsum.end(),0);
    dataType dV2sum = 0;
    vector<int>::const_iterator li = left.begin();
    vector<int>::const_iterator ri = right.begin();
    vector<bool>::const_iterator gi = groupEnd.begin();
    vector<long>::const_iterator mi = memoffset.begin();
    vector<int>::const_iterator nli = n_lookup.begin();
    dataTypeIterator vsi;
    bool isSatisfied = false;
    //
    // Iterate over all groups of neighbors, testing the balanced
    // neighborhood criterion
    //
    for (; nli != n_lookup.end(); nli++) {
      // Check to see if pixel coordinate violates group limits
      bool overEdge = false;
      for (dimIndex = 0; dimIndex < n_spatial_dims; dimIndex++,li++,ri++)
	if (pI.coord(dimIndex) < *li || pI.coord(dimIndex) >= *ri)
	  overEdge = true;
      if (overEdge)
	break;
      // Update dVsum and dV2sum with all displacements in group
      for (;; mi++,gi++) {
	//mexPrintf("(%ld) ",*mi);
	imdatatmp = imdata+*mi;
	const dataType *imdataEnd = imdatatmp+n_values;
	const dataType *imfilteredIntmp;
	for (imfilteredIntmp=imfilteredIn,vsi=dVsum.begin(); imdatatmp<imdataEnd; imdatatmp++,imfilteredIntmp++,vsi++) {
	  //mexPrintf("%g ",*imdatatmp);
	  dataType dv = *imdatatmp - *imfilteredIntmp;
	  *vsi += dv;
	  dV2sum += dv*dv;
	}
	//mexPrintf("\n");
	if (*gi) {
	  mi++;
	  gi++;
	  break;  // we are done updating this group
	}
      }
      // Test the balanced neighborhood criterion
      dVsum2 = 0;
      for (vsi = dVsum.begin(); vsi != dVsum.end(); vsi++)
	dVsum2 += *vsi * *vsi;
      isSatisfied = dVsum2 > z2*dV2sum;
      if (isSatisfied)
	break;
    }
    int n_tested = gi-groupEnd.begin();
    if (isSatisfied) {
      // The criterion was satisfied, do the backtracking
      *n = *nli;
      if (n_triggered != NULL) {
	*n_triggered = n_tested;
	n_triggered++;
      }
    } else {
      // The criterion was not satisfied, so just use the number of
      // points tested
      *n = n_tested;
      if (n_triggered != NULL) {
	*n_triggered = 0;
	n_triggered++;
      }
    }
    //mexPrintf("n = %g\n",*n);
    // If needed, compute the filtered output (the mean of all values
    // in the neighborhood)
    if (imfilteredOut != NULL) {
      //mexPrintf("Filtered:\n");
      dataType *imfilteredOutEnd = imfilteredOut+n_values;
      dataType *imfilteredOutTmp;
      fill(imfilteredOut,imfilteredOutEnd,0);
      for (int nbrI = 0; nbrI < *n; nbrI++) {
	for (imfilteredOutTmp=imfilteredOut,imdatatmp=imdata+memoffset[nbrI]; imfilteredOutTmp < imfilteredOutEnd; imfilteredOutTmp++,imdatatmp++) {
	  *imfilteredOutTmp += *imdatatmp;
	  //mexPrintf("%g ",*imdatatmp);
	}
	//mexPrintf("\n");
      }
      for ( ; imfilteredOut < imfilteredOutEnd; imfilteredOut++)
	*imfilteredOut /= *n;
    }
  }
}

// For debugging & profiling: build a stand-alone application
#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 3;
  const int n_outputs = 2;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "im",
    "memoffset",
    "coordoffset"
  };
  const char *output_names[] = {
    "n",
    "imf"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables im, memoffset, and coordoffset\n",argv[0]);
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
