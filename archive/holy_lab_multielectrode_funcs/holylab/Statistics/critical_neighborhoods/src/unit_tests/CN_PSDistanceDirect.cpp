#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerDistanceDirectC.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhoodBase.h"
#include "mat.h"
#include "MatlabIO.h"

#include <valgrind/callgrind.h>

// Compile with make

using namespace std;

int main()
{
  // Specify the types we'll be using
  typedef float DataType;
  typedef MatlabIO::MapMx<DataType>::MapMatrixType DataMatrixType;
  typedef MatlabIO::MapMx<DataType>::ColVectorType PointType;
  typedef Metric::Euclidean                        MetricType;
  typedef PointServerDistanceDirectC<DataMatrixType,MetricType> PointServerType;
  typedef AccumulatorMoments<PointType,Moments::Full> AccType;

  // Specify the p-value
  const DataType pvalue = 0.01;
  int nMin = ceil(-log(pvalue));

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("CN_PSDistanceDirect_input.mat");
  reader.load("x",X);
  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  MetricType metric;    // some metrics will need parameters passed via constructor
  PointServerType ps(X,metric);
  AccType acc(n_dims);
  CriticalNeighborhoodBase<PointServerType,AccType> cn(ps,acc);

  // Initialize outputs and auxillary storage
  Eigen::VectorXi n_checked(n_points);
  Eigen::VectorXi n(n_points);

  CALLGRIND_START_INSTRUMENTATION

  //int basePointIndex = 1246;
  int basePointIndex = 0;
  ps.basePoint(X.col(basePointIndex));
  cn.inflateAll();

  // Satisfy the critical neighborhood criterion for each point
  for (int i = 0; i < n_points; i++) {
    ps.basePoint(X.col(i));
    n[i] = cn.inflateToCriterion(pvalue,nMin);
    n_checked[i] = cn.sortOrder().size();
  }

  CALLGRIND_STOP_INSTRUMENTATION;


  // Re-run a particular point
  basePointIndex = 22;
  ps.basePoint(X.col(basePointIndex));
  int nNbrs = cn.inflateToCriterion(pvalue,nMin);
  cout << "BasePoint " << basePointIndex << ", n_checked = " << cn.sortOrder().size() << ", nNbrs = " << nNbrs << endl;
  for (int i = 0; i < cn.sortOrder().size(); i++)
    cout << ' ' << cn.sortOrder()[i];
  cout << endl;

  // Run backtracking
  cn.backtrack(nNbrs);

  // Save the output
  MatlabIO::Save writer("CN_PSDistanceDirect_output.mat");
  writer.save(n_checked,"n_checked");
  writer.save(n,"n");
  writer.saveScalar(pvalue,"pvalue");
  writer.saveScalar(basePointIndex+1,"basePointIndex");  // +1 for Matlab
  writer.save(cn.chisqNeighborhood(),"chisqNeighborhood");
  writer.save(cn.chisqPoint(),"chisqPoint");
  writer.save(cn.sortOrder(),"sortOrder");
}
