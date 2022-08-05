#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerFilter1d.h"
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
  typedef PointServerFilter1d<DataMatrixType> PointServerType;
  //typedef Metric::Euclidean                        MetricType;
  //typedef PointServerFilter1d<DataMatrixType,false,MetricType> PointServerType;

  // Specify the p-value
  const DataType pvalue = 0.01;
  int nMin = ceil(-log(pvalue));

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("PSFilter_input.mat");
  reader.load("x",X);
  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  //MetricType metric;    // some metrics will need parameters passed via constructor
  //OutlierStatistics<DataType> os(pvalue);
  PointServerType ps(X);

  int basePointIndex = 15;
  ps.basePointIndex(basePointIndex);

  for (; !ps.isAtEnd(); ps++)
    cout << "Current point index: " << ps.currentIndex() << ", current value " << ps.currentPoint() << endl;
}
