#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#else
#include "mex.h"
#endif
#include <math.h>

// Convert new and delete to using matlab memory allocation, so garbage-collection and error messages are handled properly
#include "newdelete.h"
// Code for dealing with matlab structures
#include "mat2C.cpp"
// This next file contains the routines that do all the real work
#include "msams_core.cpp"

const bool correct_offset = true;

/* msams_converge1: move point by adaptive meanshift until it
   converges, taking advantage of having a sorted list of
   neighbors. Also reports data on closest-approach to other points,
   for the purpose of assigning fates to these points.
 *
 * Syntax:
 *   nbrInfoFinal = msams_converge1(x,y)
 *   [nbrInfoFinal,status,minDist,trajectory,nbrInfoBase] = msams_converge1(x,y)
 *   [...] = msams_converge1(x,y,l2)
 *   [...] = msams_converge1(x,y,l2,lgroups)
 *   [...] = msams_converge1(...,assignment)
 *   [...] = msams_converge1(...,nbrInfoBaseIn)
 *   [...] = msams_converge1(...,options)
 * 
 * where
 *   x is a d-by-N matrix of data points in d-dimensional space (y is the "moving" probe point);
 *   y is a d-by-1 vector (the probe point) in d-dimensional space;
 *   l2 is a d-by-1 vector of squared length scales (scaling factors in the metric); if not supplied, or empty, it does not scale the metric (it is fixed at all 1s)
 *   lgroups is an optional d-by-1 vector of integers between 0 and d. All coordinates of l2 with the same integer value will be "tied together," meaning they will have the same value (which, after updates, is set to the mean of the individual values). In other words, lgroups = [1 2 1 2 2] specifies that there will be only 2 unique values for the coordinates of l2, one value for the first and third coordinates, and the other value for the second, fourth, and fifth coordinates.
 *   assignment (if provided) is a 1-by-N vector of integers, one for each point in x; 0 means unassigned, any positive integer is the "assignment index" for the corresponding point in x. If y's nearest neighbor is assigned (assignment > 0), the algorithm will truncate early, returning the assignment in the "status" output.
 *   options is a structure which may contain the following fields
 *     factor (default 3)
 *     min_to_check (default d+1)
 *     backtrack (default true)
 *     n_min (default 0)
 *     n_threads (default # of cores)
 *
 * On output:
 *   yfinal is the final position of this probe point
 *   l2final is the mean square size of the neighborhood (=final set of scaling factors)
 *   n is the # points in yfinal's neighborhood
 *   nbrWeight is a 1-by-N matrix with entries lying between 0 and 1
 *     inclusive, indicating the maximum "degree of closeness" that each
 *     point in x had to the probe point as it moved towards its peak
 *     (0 = never a member of the neighborhood, 1 = right on top of
 *       the probe point)
 *   status is a structure with the following fields:
 *     converged: false if didn't converge in the # of iterations used
 *     assignment: index of the nearest-neightbor's assignment (if nonzero, converged is likely to be false, because it exits early)
 *     n_iter: a two-vector containing the number of "full" steps (in the first position) and the number of neighborhood-based steps (in the second position)
 *   nbrInfoOut: the balanced_neighborhood structure. Note this is NOT the sorting of neighbors relative to the returned position, it's the last "direct" (full) ordering of the neighbors.
 */

// Utility functions
template <class T> 
void matOptionalScalarField2C(const mxArray *m,const char *name,T *v);
void mat2C(const mxArray *mxOptions,optionStruct &ops);
template <class dataType>
    void mat2C(const mxArray *nbrInfo,balanced_neighborhood<dataType> &nbr);

const char *status_fields[] = {
  "converged",
  "assignment",
  "n_iter"};
int n_status_fields = 3;

const char *balanced_neighborhood_fields[] = {
  "ybase",
  "l2",
  "neighborLabel",
  "neighborDist",
  "n_exceed",
  "n",
  "nbrhoodMean",
  "nbrhoodMSDisp",
  "RMSDist"};
  int n_balanced_neighborhood_fields = 9;

const char *trajectory_fields[] = {
  "y",
  "l2",
  "n",
  "nearestNeighbor"};
int n_trajectory_fields = 4;

// This is the "inner" matlab wrapper. It is called by mexFunction
// (below), but this one does all the real work.  It's templated so
// that it can work with a variety of data types.
template <class dataType,class intType>
void msams_converge1_wrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[],mxClassID mat_dataType,mxClassID mat_intType)
{
  data_points<dataType> points;
  optionStruct ops;
  balanced_neighborhood<dataType> nbrbase;
  balanced_neighborhood<dataType> nbr;
  const mxArray *curarg;
  mxArray *curfield;
  int argindx,curfieldN;
  int dimIterator,ptIterator;

  // Parse the inputs
  // Get the data points
  curarg = prhs[0];
  // We don't have to check the type of the data points, because
  // that's done in mexFunction
  points.x = (dataType*) mxGetData(curarg);
  points.d = mxGetM(curarg);
  points.N = mxGetN(curarg);

  // Get y
  curarg = prhs[1];
  if (mxGetClassID(curarg) != mat_dataType || mxIsComplex(curarg))
    mexErrMsgTxt("msams_converge1: y must be of the same type as the data");
  if (mxGetNumberOfElements(curarg) != points.d)
    mexErrMsgTxt("msams_converge1: y must have d elements");
  dataType *y = (dataType*) mxGetData(curarg);

  // Parse the remaining arguments. Since the remaining inputs are
  // optional, we have to deduce which is which
  ops.min_to_check = points.d+1;
  dataType *l2 = NULL;
  intType *lgroups = NULL;
  intType *assignment = NULL;
  const mxArray *optionsP = NULL;
  nbrbase.allocate(points.d,points.N);
  for (argindx = 2; argindx < nrhs; argindx++) {
    curarg = prhs[argindx];
    if (mxGetClassID(curarg) == mat_dataType && mxGetM(curarg) == points.d && mxGetN(curarg) == 1 && l2 == NULL) {
      // It's "l2"
      if (mxIsComplex(curarg))
	mexErrMsgTxt("msams_converge1: l2 must be a column vector of the same type as the data");
      l2 = (dataType*) mxGetData(curarg);
    }
    else if (mxGetClassID(curarg) == mat_intType && mxGetM(curarg) == points.d && mxGetN(curarg) == 1 && l2 != NULL && lgroups == NULL) {
      // It's "lgroups"
      if (mxIsComplex(curarg))
        mexErrMsgTxt("msams_converge1: lgroups must be a column vector of the same type as the data");
      lgroups = (intType*) mxGetData(curarg);
    }
    else if (mxGetClassID(curarg) == mat_intType && !mxIsComplex(curarg) && mxGetM(curarg) == 1 && mxGetN(curarg) == points.N && assignment == NULL) {
      // It's "assignment"
      assignment = (intType*) mxGetData(curarg);
    } 
    else if (mxIsStruct(curarg)) {
      curfieldN = mxGetFieldNumber(curarg,"neighborLabel");
      if (curfieldN == -1) {
	// It's "options"
	mat2C(curarg,ops);
        optionsP = curarg;
      }
      else {
	// It's "nbrInfoBaseIn"
	mat2C(curarg,nbrbase);
        if (correct_offset)
          for (ptIterator = 0; ptIterator < nbr.n_valid; ptIterator++)
            nbrbase.neighborLabel[ptIterator]--;  // correct for matlab unit offset
      }
    }
    else if (mxIsEmpty(curarg))
      continue;
    else
      mexErrMsgIdAndTxt("msams:converge","Error parsing argument %d",argindx+1);
  }

  // Check options for sanity
  if (ops.min_to_check > points.N)
    mexWarnMsgIdAndTxt("msams:options","options.min_to_check is larger than the number of points. This may be due to high-dimensionality. It is recommended that you manually set options.min_to_check.");

  if (ops.min_to_check < ops.n_min)
    ops.min_to_check = ops.n_min;
  ops.min_to_check = (points.N >= ops.min_to_check) ? ops.min_to_check : points.N;  // set to min(min_to_check,N)
  if (ops.n_min > ops.min_to_check)
    ops.n_min = ops.min_to_check;
  //ops.n_threads = (ops.n_cpus < ops.n_threads) ? ops.n_cpus : ops.n_threads;

  // If l2 is empty, set up a default
  bool l2empty = (l2 == NULL);
  if (l2empty) {
    l2 = (dataType*) mxMalloc(points.d*sizeof(dataType));
    for (dimIterator = 0; dimIterator < points.d; dimIterator++)
      l2[dimIterator] = 1.0;
    ops.variable_metric = false;
  }
  // Set up default for lgroups
  bool free_lgroups = false;
  if (ops.variable_metric && lgroups == NULL) {
    lgroups = new intType[points.d];
    free_lgroups = true;
    for (dimIterator = 0; dimIterator < points.d; dimIterator++)
      lgroups[dimIterator] = dimIterator;
  }
  // Set up the l2 consolidator
  l2_operations<dataType,intType> l2func;
  l2func.initialize(lgroups,points.d);
  if (optionsP != NULL && !l2empty) {
    curfield = mxGetField(optionsP,0,"l2min");
    if (curfield != NULL) {
      if (mxGetClassID(curfield) != mat_dataType || mxGetM(curfield) != points.d || mxGetN(curfield) != 1)
        mexErrMsgTxt("options.l2min must be of the same type and shape as a data point");
      l2func.set_minimum((dataType*) mxGetData(curfield));
    }
  }

  // Allocate the outputs or create storage for the algorithm
  // status
  mxArray *statusP = NULL;
  int *niterP = NULL;
  if (nlhs > 1) {
    statusP = mxCreateStructMatrix(1,1,n_status_fields,status_fields);
    plhs[1] = statusP;
    curfield = mxCreateNumericMatrix(1,2,mxINT32_CLASS,mxREAL);
    mxSetField(statusP,0,"n_iter",curfield);
    niterP = (int *) mxGetData(curfield);
  }
  // minDist
  dataType *minDist = NULL;
  if (nlhs > 2) {
    plhs[2] = mxCreateNumericMatrix(1,points.N,mat_dataType,mxREAL);
    minDist = (dataType*) mxGetData(plhs[2]);
  }
  // Set up the trajectory storage structure
  trajectory_operations<dataType> trajectory;
  mxArray *trajP = NULL;
  if (nlhs > 3) {
    trajP = matlab_allocate_trajectory(trajectory,points.d,ops.iter_max,mat_dataType);
    plhs[3] = trajP;
  }
  else
    trajectory.allocate(points.d,ops.iter_max);
  nbr.allocate(points.d,points.N);

  // Call the main algorithm
  bool converged;
  bool is_assigned;
  try {
    converged = msams_converge_core(points,nbrbase,y,l2,l2func,assignment,nbr,minDist,is_assigned,niterP,trajectory,ops);
  }
  catch (std::runtime_error &e) {
      mexErrMsgTxt(e.what());
  }

  // Pack the outputs
  int nearestNbrIndex = nbr.neighborLabel[0];
  // nbrInfo
  mxArray *nbrInfoP = matlab_allocate_bn(points.d,nbr.n_valid,mat_dataType);
  plhs[0] = nbrInfoP;
  // offset neighborLabel by 1 (matlab is unit-offset)
  if (correct_offset)
    for (ptIterator = 0; ptIterator < nbr.n_valid; ptIterator++)
      nbr.neighborLabel[ptIterator]++;
  C2mat(nbr,nbrInfoP);
  //status
  if (statusP != NULL) {
    curfield = mxCreateLogicalScalar(converged);
    mxSetField(statusP,0,"converged",curfield);
    if (assignment != NULL && is_assigned)
      curfield = mxCreateDoubleScalar(double(assignment[nearestNbrIndex]));
    else
      curfield = mxCreateDoubleScalar(0);
    mxSetField(statusP,0,"assignment",curfield);
  }
  // minDist is already handled
  // trajectory
  if (trajP != NULL) {
    // offset nearestNeighbor by 1 (for matlab unit-offset)
    for (ptIterator = 0; ptIterator < trajectory.traj_length; ptIterator++)
      trajectory.nntraj[ptIterator]++;
    matlab_reallocate_trajectory(trajP,trajectory.traj_length);
  }
  // nbrInfoBaseOut
  if (nlhs > 4) {
    mxArray *nbrInfoBaseP = matlab_allocate_bn(points.d,points.N,mat_dataType);
    plhs[4] = nbrInfoBaseP;
    if (correct_offset)
      for (ptIterator = 0; ptIterator < nbr.n_valid; ptIterator++)
        nbrbase.neighborLabel[ptIterator]++;
    C2mat(nbrbase,nbrInfoBaseP);
  }

   // Free temporary storage
  if (l2empty)
    mxFree(l2);
  if (free_lgroups)
    mxFree(lgroups);

  return;
}

// The "outer" matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nrhs < 2 || nrhs > 7)
    mexErrMsgTxt("msams_converge1: requires 2-7 inputs");
  if (nlhs < 1 || nlhs > 5)
    mexErrMsgTxt("msams_converge1: requires 1-5 outputs");

  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfDimensions(curarg) != 2)
    mexErrMsgTxt("msams_converge1: x must be a real matrix");

  if (mxIsDouble(curarg))
    msams_converge1_wrapper<double,double>(nlhs,plhs,nrhs,prhs,mxDOUBLE_CLASS,mxDOUBLE_CLASS);
  else if (mxIsSingle(curarg))
    msams_converge1_wrapper<float,double>(nlhs,plhs,nrhs,prhs,mxSINGLE_CLASS,mxDOUBLE_CLASS);
  else
    mexErrMsgTxt("msams_converge1: x must be a single- or double-precision");

  return;
}

void mat2C(const mxArray *mxOptions,optionStruct &ops)
{
  matOptionalScalarField2C(mxOptions,"factor",&(ops.factor));
  matOptionalScalarField2C(mxOptions,"min_to_check",&(ops.min_to_check));
  matOptionalScalarField2C(mxOptions,"backtrack",&(ops.backtrack));
  matOptionalScalarField2C(mxOptions,"n_min",&(ops.n_min));
  matOptionalScalarField2C(mxOptions,"iter_max",&(ops.iter_max));
  matOptionalScalarField2C(mxOptions,"variable_metric",&(ops.variable_metric));
  matOptionalScalarField2C(mxOptions,"tol",&(ops.tol));
  matOptionalScalarField2C(mxOptions,"n_threads",&(ops.n_threads));
}

template <class dataType>
    void mat2C(const mxArray *nbrInfo,balanced_neighborhood<dataType> &nbr)
{
  // nbr must already be allocated, and the data are copied-by-value into nbr
  int n;
  matField2Carray(nbrInfo,"ybase",nbr.ybase,nbr.d,nbr.d);
  matField2Carray(nbrInfo,"l2",nbr.l2,nbr.d,nbr.d);
  n = matField2Carray(nbrInfo,"neighborLabel",nbr.neighborLabel,nbr.N,0);
  matField2Carray(nbrInfo,"neighborDist",nbr.neighborDist,n,n);
  nbr.n_valid = n;
  matScalarField2C(nbrInfo,"n_exceed",&nbr.n_exceed);
  matScalarField2C(nbrInfo,"n",&nbr.n);
  matField2Carray(nbrInfo,"nbrhoodMean",nbr.nbrhoodMean,nbr.d,nbr.d);
  matField2Carray(nbrInfo,"nbrhoodMSDisp",nbr.nbrhoodMSDisp,nbr.d,nbr.d);
  matScalarField2C(nbrInfo,"RMSDist",&nbr.RMSDist);
}

mxArray *matlab_allocate_bn(int d,int n_valid,mxClassID mat_dataType)
{
  mxArray *curfield;

  mxArray *nbrInfoP = mxCreateStructMatrix(1,1,n_balanced_neighborhood_fields,balanced_neighborhood_fields);
  curfield = mxCreateNumericMatrix(d,1,mat_dataType,mxREAL);
  mxSetField(nbrInfoP,0,"ybase",curfield);
  curfield = mxCreateNumericMatrix(d,1,mat_dataType,mxREAL);
  mxSetField(nbrInfoP,0,"l2",curfield);
  curfield = mxCreateNumericMatrix(1,n_valid,mxINT32_CLASS,mxREAL);
  mxSetField(nbrInfoP,0,"neighborLabel",curfield);
  curfield = mxCreateNumericMatrix(1,n_valid,mat_dataType,mxREAL);
  mxSetField(nbrInfoP,0,"neighborDist",curfield);
  curfield = mxCreateDoubleScalar(-1);
  mxSetField(nbrInfoP,0,"n_exceed",curfield);
  curfield = mxCreateDoubleScalar(-1);
  mxSetField(nbrInfoP,0,"n",curfield);
  curfield = mxCreateDoubleScalar(-1);
  mxSetField(nbrInfoP,0,"RMSDist",curfield);
  curfield = mxCreateNumericMatrix(d,1,mat_dataType,mxREAL);
  mxSetField(nbrInfoP,0,"nbrhoodMean",curfield);
  curfield = mxCreateNumericMatrix(d,1,mat_dataType,mxREAL);
  mxSetField(nbrInfoP,0,"nbrhoodMSDisp",curfield);

  return nbrInfoP;
}

template <class dataType>
    mxArray *matlab_allocate_trajectory(trajectory_operations<dataType> &trajectory,int d,int iter_max,mxClassID mat_dataType)
{
  mxArray *curfield;
  mxArray *trajP = mxCreateStructMatrix(1,1,n_trajectory_fields,trajectory_fields);
  curfield = mxCreateNumericMatrix(d,iter_max,mat_dataType,mxREAL);
  mxSetField(trajP,0,"y",curfield);
  trajectory.ytraj = (dataType*) mxGetData(curfield);
  curfield = mxCreateNumericMatrix(d,iter_max,mat_dataType,mxREAL);
  mxSetField(trajP,0,"l2",curfield);
  trajectory.l2traj = (dataType*) mxGetData(curfield);
  curfield = mxCreateNumericMatrix(1,iter_max,mxINT32_CLASS,mxREAL);
  mxSetField(trajP,0,"n",curfield);
  trajectory.ntraj = (int32_t*) mxGetData(curfield);
  curfield = mxCreateNumericMatrix(1,iter_max,mxINT32_CLASS,mxREAL);
  mxSetField(trajP,0,"nearestNeighbor",curfield);
  trajectory.nntraj = (int32_t*) mxGetData(curfield);
  trajectory.d = d;

  return trajP;
}

template <class dataType>
    void C2mat(balanced_neighborhood<dataType> &nbr,mxArray *nbrInfoP)
{
  Carray2matField(nbrInfoP,"ybase",nbr.ybase,nbr.d);
  Carray2matField(nbrInfoP,"l2",nbr.l2,nbr.d);
  Carray2matField(nbrInfoP,"neighborLabel",nbr.neighborLabel,nbr.n_valid);
  Carray2matField(nbrInfoP,"neighborDist",nbr.neighborDist,nbr.n_valid);
  C2matScalarField(nbrInfoP,"n_exceed",(double) nbr.n_exceed);
  C2matScalarField(nbrInfoP,"n",(double) nbr.n);
  C2matScalarField(nbrInfoP,"RMSDist",(double) nbr.RMSDist);
  Carray2matField(nbrInfoP,"nbrhoodMean",nbr.nbrhoodMean,nbr.d);
  Carray2matField(nbrInfoP,"nbrhoodMSDisp",nbr.nbrhoodMSDisp,nbr.d);
}

void matlab_reallocate_trajectory(mxArray *trajP,int iter_max)
{
  mxArray *curfield;
  curfield = mxGetField(trajP,0,"y");
  mxSetN(curfield,iter_max);
  curfield = mxGetField(trajP,0,"l2");
  mxSetN(curfield,iter_max);
  curfield = mxGetField(trajP,0,"n");
  mxSetN(curfield,iter_max);
  curfield = mxGetField(trajP,0,"nearestNeighbor");
  mxSetN(curfield,iter_max);
}

// For debugging & profiling: build a stand-alone application
#ifdef MAIN
//[yfinal,l2final,n,nbrWeight] = msams_converge1(x,y,l2,options)
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 6;
  const int n_outputs = 4;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "x",
    "y",
    "l2",
    "lgroups",
    "assignment",
    "options"};
  const char *output_names[] = {
    "nbrInfo",
    "status",
    "minDist",
    "trajectory"};

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables x, y, l2, assignment, and options\n",argv[0]);
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
  //if (debug)
  //  printf("About to save\n");
  mat_save_variables(fileout,output_names,n_outputs,output);
  //if (debug)
  printf("All done saving\n");

  return EXIT_SUCCESS;
}
#endif
