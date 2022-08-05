#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerDistanceDirectC.h"
#include "mat.h"
#include "MatlabIO.h"
#include "OutlierStatistics.h"

#include <valgrind/callgrind.h>

// Compile with make

using namespace std;

int main()
{
  // Specify the types we'll be using
  typedef double DataType;
  typedef MatlabIO::MapMx<DataType>::MapMatrixType DataMatrixType;
  typedef MatlabIO::MapMx<DataType>::ColVectorType PointType;
  typedef Metric::Euclidean                        MetricType;
  typedef PointServerDistanceDirectC<DataMatrixType,MetricType> PointServerType;

  // Specify the p-value
  const DataType pvalue = 0.01;
  int nMin = ceil(-log(pvalue));

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("PSDistanceDirect_input.mat");
  reader.load("x",X);
  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  MetricType metric;    // some metrics will need parameters passed via constructor
  OutlierStatistics<DataType> os(pvalue);
  PointServerType ps(X,metric,os);

  // Initialize outputs and auxillary storage
  Eigen::VectorXi n_checked(n_points);
  Eigen::VectorXi n(n_points);

  int basePointIndex = 0;
  ps.basePoint(X.col(basePointIndex));

  cout << "The total number of points in the data set is " << ps.Ntot()
       << ",\nand the number before the first outlier is " << ps.N() << endl;

  /*
  // Save the output
  MatlabIO::Save writer("CN_PSDistanceDirect_output.mat");
  writer.save(n_checked,"n_checked");
  writer.save(n,"n");
  writer.saveScalar(pvalue,"pvalue");
  writer.saveScalar(basePointIndex+1,"basePointIndex");  // +1 for Matlab
  writer.save(cn.chisqNeighborhood(),"chisqNeighborhood");
  writer.save(cn.chisqPoint(),"chisqPoint");
  writer.save(cn.sortOrder(),"sortOrder");
  */
}
