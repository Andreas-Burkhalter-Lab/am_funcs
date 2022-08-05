#include "mex.h"
#include "PointServerFilter1d.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhoodBase.h"
#include "MatlabIO.h"

using namespace Eigen;

// See cn_filter.m for the syntax

// This performs the actual filtering computation. It's templated so that we can call it for the whole zoo of possibilities.
template <typename PSType,typename AccType>
void run_filter(PSType& ps,AccType& acc,NeighborhoodStatistics<typename AccType::Scalar>& nstats,NeighborhoodHistory& history,typename PSType::MatrixType& xf,double* pn,double* pniter) {
  CriticalNeighborhoodBase<PSType,AccType> cn(ps,acc,nstats);
  for (int i = 0; i < ps.Ntot(); i++) {
    ps.basePointIndex(i);
    CriticalNeighborhoodAlgorithms::flowToPeak(cn,history);
    cn.backtrackPeak();
    xf.col(i) = acc.mean();
    if (pn != NULL)
      pn[i] = acc.n();
    if (pniter != NULL)
      pniter[i] = history.n_history().size();
  }
}

// This is the "inner" matlab wrapper. It does all the real argument parsing.
// It's templated so that it can work with a variety of data types.
template <class Scalar>
void mexWrapper(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
  typedef Matrix<Scalar,Dynamic,1> PointType;
  typedef Map<MatrixType> MapType;
  // Types that represent the mix-and-match possibilities
  typedef PointServerFilter1d<MapType,true> PSAcausalType;
  typedef PointServerFilter1d<MapType,false> PSCausalType;

  int d,N;
  Scalar pvalue;
  bool isAcausal;
  Moments::CovarianceModel covarianceModel = Moments::Isotropic;

  // Get the values
  MapType x(NULL,0,0);
  MatlabIO::mx2eigen(prhs[0],x);
  d = x.rows();
  N = x.cols();
  if (N == 1)
    mexErrMsgTxt("Cannot filter a signal that is only a single time point. Did you mean to supply the transpose?");

  // Get the p-value
  if (mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgTxt("pvalue must be a scalar");
  pvalue = mxGetScalar(prhs[1]);

  // Get the causality flag
  isAcausal = false;
  int curarg = 2;
  if (nrhs > curarg) {
    if (mxGetNumberOfElements(prhs[curarg]) != 1)
      mexErrMsgTxt("isAcausal must be a scalar");
    isAcausal = (bool) mxGetScalar(prhs[curarg]);
  }

  // Get the covariance model
  curarg = 3;
  if (nrhs > curarg) {
    if (!mxIsChar(prhs[curarg]))
      mexErrMsgTxt("covarianceModel must be a string");
    mxChar *pc = mxGetChars(prhs[curarg]);
    if (*pc == 'i')
      ;  // Default value, do nothing
    else if (*pc == 'd')
      covarianceModel = Moments::Diagonal;
    else if (*pc == 'f')
      covarianceModel = Moments::Full;
    else
      mexErrMsgIdAndTxt("cn:covariance","covarianceModel %s not recognized",pc);
  }

  // Create the output storage
  // xf
  plhs[0] = MatlabTraits<Scalar>::mxAllocator(2,mxGetDimensions(prhs[0]));
  MapType xf(NULL,0,0);
  MatlabIO::mx2eigen(plhs[0],xf);

  // n
  double *pn = NULL;
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
    pn = mxGetPr(plhs[1]);
  }

  // n_iter
  double *pniter = NULL;
  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
    pniter = mxGetPr(plhs[2]);
  }

  // Allocate the storage for neighborhood statistics & history
  NeighborhoodStatistics<Scalar> nstats(pvalue,0,covarianceModel);
  NeighborhoodHistory history;


  // Run the algorithm
  if (covarianceModel == Moments::Isotropic) {
    typedef AccumulatorMoments<PointType,Moments::Isotropic> AccType;
    // Isotropic case
    AccType acc(d);
    if (isAcausal) {
      typedef PSAcausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    } else {
      typedef PSCausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    }
  } else if (covarianceModel == Moments::Diagonal) {
    // Diagonal case
    typedef AccumulatorMoments<PointType,Moments::Diagonal> AccType;
    AccType acc(d);
    if (isAcausal) {
      typedef PSAcausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    } else {
      typedef PSCausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    }
  } else {
      // Full case
    typedef AccumulatorMoments<PointType,Moments::Full> AccType;
    AccType acc(d);
    if (isAcausal) {
      typedef PSAcausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    } else {
      typedef PSCausalType PSType;
      PSType ps(x);
      run_filter(ps,acc,nstats,history,xf,pn,pniter);
    }
  }
}


// The overall matlab wrapper. All this does is act as a switchyard
// for choosing the single-precision or double-precision templated code.
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;

  if (nlhs == 0)
    return;
  else if (nlhs > 3)
    mexErrMsgTxt("Only one to three output arguments are provided");
  if (nrhs < 2 || nrhs > 4)
    mexErrMsgTxt("Requires 2-4 inputs: x, pvalue, iscausal (optional), covarianceModel (optional)");

  // Determine the data type from the first argument
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("x must be a real numeric array");

  if (mxIsDouble(curarg))
    mexWrapper<double>(nlhs,plhs,nrhs,prhs);
  else if (mxIsSingle(curarg))
    mexWrapper<float>(nlhs,plhs,nrhs,prhs);
  else
    mexErrMsgTxt("x must be a single- or double-precision");
    
  return;
}


