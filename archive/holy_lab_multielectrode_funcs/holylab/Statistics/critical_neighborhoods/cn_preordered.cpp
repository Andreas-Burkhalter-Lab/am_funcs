#include "mex.h"
#include "PointServerPreorderedC.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhoodAlgorithms.h"
#include "MatlabIO.h"

using namespace Eigen;

const int matlab_offset = 1;

// See cn_preordered.m for the syntax

// This is the "inner" matlab wrapper. It does all the real work.
// It's templated so that it can work with a variety of data types.

template <class Scalar>
void mexWrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
  typedef Matrix<Scalar,Dynamic,1> PointType;
  typedef Map<MatrixType> MapType;
  typedef double OrderScalar;
  typedef std::vector<OrderScalar> OrderType;

  int d,covariance_model;
  Scalar scalartmp;
  std::vector<Scalar> T2;
  OrderType sortOrderIn;
  std::vector<bool> tiesPrevious;
  std::vector<int> sortOrder;
  bool useWeight;

  // Parse the inputs
  // Get the displacements
  MapType dX(NULL,0,0);
  MatlabIO::mx2eigen(prhs[0],dX);
  d = dX.rows();

  // Get the weights
  // Determine whether weights are supplied by whether the 3rd argument
  // is a string
  useWeight = (nrhs > 2 && !mxIsChar(prhs[2]));
  MapType weight(NULL,0,0);
  if (useWeight) {
    MatlabIO::mx2eigen(prhs[1],weight);
    if (weight.cols() != dX.cols())
      mexErrMsgTxt("Weight must have the same number of columns as the data matrix");
    if (weight.rows() != 1)
      mexErrMsgTxt("Weight must have a single row");
  }

  // Get the sortOrder
  MatlabIO::mx2vector(prhs[1+useWeight],sortOrderIn);
  //  if (sortOrderIn.size() != dX.cols())
  //    mexErrMsgTxt("sortOrder must have the same number of elements as there are points");

  // Get the covariance model
  covariance_model = 0;
  int curarg = 2+useWeight;
  if (nrhs > curarg) {
    if (!mxIsChar(prhs[curarg]))
      mexErrMsgTxt("covarianceModel must be a string");
    mxChar *pc = mxGetChars(prhs[curarg]);
    if (*pc == 'i')
      covariance_model = 0;
    else if (*pc == 'd')
      covariance_model = 1;
    else if (*pc == 'f')
      covariance_model = 2;
    else
      mexErrMsgIdAndTxt("cn:covariance","covarianceModel %s not recognized",pc);
  }
    
  // Create the base point (all zeros)
  PointType bp(d);
  bp.setZero();

  if (covariance_model == 0) {
    // Isotropic case
    AccumulatorMoments<PointType,Moments::Isotropic> acc(d);
    if (useWeight) {
      PointServerPreorderedCWeighted<MapType,MapType,OrderType,1> ps(dX,weight,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    } else {
      PointServerPreorderedC<MapType,OrderType,1> ps(dX,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    }
  } else if (covariance_model == 1) {
    // Diagonal case
    AccumulatorMoments<PointType,Moments::Diagonal> acc(d);
    if (useWeight) {
      PointServerPreorderedCWeighted<MapType,MapType,OrderType,1> ps(dX,weight,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    } else {
      PointServerPreorderedC<MapType,OrderType,1> ps(dX,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    }
  } else {
      // Full case
    AccumulatorMoments<PointType,Moments::Full> acc(d);
    if (useWeight) {
      PointServerPreorderedCWeighted<MapType,MapType,OrderType,1> ps(dX,weight,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    } else {
      PointServerPreorderedC<MapType,OrderType,1> ps(dX,sortOrderIn);
      ps.basePoint(bp);
      CriticalNeighborhoodAlgorithms::inflateAll(T2,sortOrder,tiesPrevious,acc,ps);
    }
  }

  plhs[0] = MatlabIO::vector2mx(T2);
}


// The overall matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nlhs == 0)
    return;
  else if (nlhs > 2)
    mexErrMsgTxt("Only one or two output arguments are provided");
  if (nrhs < 2 || nrhs > 4)
    mexErrMsgTxt("Requires 2-4 inputs: dX, weight (optional), sortOrder, covariance_model (optional)");

  // Determine the data type from the first argument
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("dX must be a real numeric array");

  if (mxIsDouble(curarg))
    mexWrapper<double>(nlhs,plhs,nrhs,prhs);
  else if (mxIsSingle(curarg))
    mexWrapper<float>(nlhs,plhs,nrhs,prhs);
  else
    mexErrMsgTxt("dX must be a single- or double-precision");
    
  return;
}


