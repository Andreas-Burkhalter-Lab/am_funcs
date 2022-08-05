#include <iostream>
#include <vector>
#include "PointServerDistanceDirectC.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhoodBase.h"
#include "mat.h"
#include "MatlabIO.h"

#include <valgrind/callgrind.h>

// Compile with make

using namespace std;

template <typename PointServer,typename Accumulator>
void clusterAllFlow(vector<int>& mapsTo,vector<int>& n,PointServer& ps,Accumulator& acc,double pvalue)
{
  int n_points = ps.N();
  int nNbrs;
  int nMin = ceil(-log(pvalue));
  NeighborhoodHistory history;
  int i,j;
  CriticalNeighborhoodBase<PointServer,Accumulator> cn(ps,acc);
  CriticalNeighborhoodAlgorithms::FlowResults result;
  vector<bool> isAtMax(n_points);

  // Move each probe point until it maps to another, or until it reaches a peak
  cout << "First pass...";
  for (i = 0; i < n_points; i++) {  // loop over probe points
    ps.basePointIndex(i);
    result = CriticalNeighborhoodAlgorithms::flowToNeighbor(cn,history,pvalue,nMin);
    mapsTo[i] = cn.sortOrder()[0];
    n[i] = result.nNbrs;
    isAtMax[i] = result.isAtMax;
  }
  cout << "done." << endl;
  for (i = 0; i < n_points; i++)
    cout << i << '/' << mapsTo[i] << "(" << n[i] << ")  ";
  cout << endl;

  // Flow the map
  vector<int> mapsTo0 = mapsTo;   // preserve the original
  CriticalNeighborhoodAlgorithms::flowMap(mapsTo,n);
  vector<int> mapsToOld;

  // Reflow until we get consistency
  bool isChanged = true;
  while (isChanged) {
    mapsToOld = mapsTo;
    // Check points that appear to map to themselves, by flowing them all the
    // way to their peak
    for (i = 0; i < n_points; i++) {
      if (mapsTo[i] == i && !isAtMax[i]) {
	cout << i;
	ps.basePointIndex(i);
	nNbrs = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue,nMin);
	mapsTo0[i] = cn.sortOrder()[0];
	n[i] = nNbrs;
	isAtMax[i] = true;
	cout << " -> " << cn.sortOrder()[0] << " (" << nNbrs << "), " << acc.mean().transpose() << endl;
      }
    }
    mapsTo = mapsTo0;
    CriticalNeighborhoodAlgorithms::flowMap(mapsTo,n);
    for (i = 0; i < n_points; i++)
      if (mapsTo[i] != mapsToOld[i])
	break;
    if (i == n_points)
      isChanged = false;
  }
}


int main()
{
  // Specify the types we'll be using
  typedef double DataType;
  typedef MatlabIO::MapMx<DataType>::MapMatrixType DataMatrixType;
  typedef MatlabIO::MapMx<DataType>::ColVectorType PointType;
  typedef Metric::Euclidean                        MetricType;
  typedef PointServerDistanceDirectC<DataMatrixType,MetricType> PointServerType;

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("cluster_all_input.mat");
  reader.load("x",X);
  // Load the parameters that affect the outcome
  mxArray* pCovModel;
  pCovModel = reader.load("covariance_model");
  int covModel = (int) mxGetScalar(pCovModel);
  mxArray* pPValue;
  pPValue = reader.load("pvalue");
  double pvalue = mxGetScalar(pPValue);

  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  MetricType metric;    // some metrics will need parameters passed via constructor
  PointServerType ps(X,metric);

  // Initialize outputs and auxillary storage
  vector<int> mapsTo(n_points);
  vector<int> n(n_points);

  if (covModel == 0) {
    AccumulatorMoments<PointType,Moments::Isotropic> acc(n_dims);
    clusterAllFlow(mapsTo,n,ps,acc,pvalue);
  } else if (covModel == 1) {
    //Eigen::Matrix<int,PointType::RowsAtCompileTime,1> group(n_dims);
    //group.setConstant(0);
    //AccumulatorMoments<PointType,Moments::Diagonal> acc(group);
    AccumulatorMoments<PointType,Moments::Diagonal> acc(n_dims);
    clusterAllFlow(mapsTo,n,ps,acc,pvalue);
  } else {
    AccumulatorMoments<PointType,Moments::Full> acc(n_dims);
    clusterAllFlow(mapsTo,n,ps,acc,pvalue);
  }

  // Save the output
  MatlabIO::Save writer("cluster_all_output.mat");
  writer.save(n,"n");
  writer.save(mapsTo,"map");
}
