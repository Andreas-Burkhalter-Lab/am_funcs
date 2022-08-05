#include <iostream>
#include <vector>
#include "PointServerDistanceDirectC.h"
#include "AccumulatorMoments.h"
#include "CriticalNeighborhood.h"
#include "mat.h"
#include "MatlabIO.h"

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

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("cluster_all_input.mat");
  reader.load("x",X);
  // Load the parameters that affect the outcome
  mxArray* pPValue;
  pPValue = reader.load("pvalue");
  double pvalue = mxGetScalar(pPValue);

  int n_dims = X.rows();
  int n_points = X.cols();
  int i;

  // Initialize the variables needed for critical neighborhood analysis
  MetricType metric;    // some metrics will need parameters passed via constructor
  PointServerType ps(X,metric);
  vector<DataType> chisqN,chisqP;
  vector<int> sortOrder;

  ps.basePointIndex(0);

  // Compare the different accumulators
  AccumulatorMoments<PointType,Moments::Isotropic> accIso(n_dims);
  CriticalNeighborhood::chisqNeighbors(chisqN,chisqP,sortOrder,accIso,ps);
  for (i = 0; i < n_points; i++)
    cout << i << " (" << chisqN[i] << '/' << chisqP[i] << ")   ";
  cout << endl;

  Eigen::Matrix<int,PointType::RowsAtCompileTime,1> group(n_dims);
  group.setConstant(0);
  AccumulatorMoments<PointType,Moments::Diagonal> accDiagGroup(group);
  CriticalNeighborhood::chisqNeighbors(chisqN,chisqP,sortOrder,accDiagGroup,ps);
  for (i = 0; i < n_points; i++)
    cout << i << " (" << chisqN[i] << '/' << chisqP[i] << ")   ";
  cout << endl;

  AccumulatorMoments<PointType,Moments::Diagonal> accDiag(n_dims);
  CriticalNeighborhood::chisqNeighbors(chisqN,chisqP,sortOrder,accDiag,ps);
  for (i = 0; i < n_points; i++)
    cout << i << " (" << chisqN[i] << '/' << chisqP[i] << ")   ";
  cout << endl;

  AccumulatorMoments<PointType,Moments::Full> accFull(n_dims);
  CriticalNeighborhood::chisqNeighbors(chisqN,chisqP,sortOrder,accFull,ps);
  for (i = 0; i < n_points; i++)
    cout << i << " (" << chisqN[i] << '/' << chisqP[i] << ")   ";
  cout << endl;

}
