#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerPreorderedC.h"
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
  typedef MatlabIO::MapMx<DataType>::MapMatrixType    DataMatrixType;
  typedef MatlabIO::MapMx<DataType>::ColVectorType    PointType;
  typedef MatlabIO::MapMx<DataType>::MapColVectorType MxPointType;
  typedef Eigen::Map<Eigen::RowVectorXi>              OrderType;
  typedef PointServerPreorderedC<DataMatrixType,OrderType,1> PointServerType;
  typedef AccumulatorMoments<PointType,Moments::Isotropic> AccType;

  // Specify the p-value
  const DataType pvalue = 0.001;

  // Load the data
  DataMatrixType X(NULL,0,0);
  OrderType sortOrder(NULL,0);
  MxPointType bp(NULL,0);
  MatlabIO::Load reader("CN_PSPreordered_input.mat");
  reader.load("x",X);
  reader.load("sortOrder",sortOrder);
  reader.load("y",bp);

  cout << "read sortOrder size " << sortOrder.size() << endl;

  CALLGRIND_START_INSTRUMENTATION

  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  PointServerType ps(X,sortOrder);
  AccType acc(n_dims);
  CriticalNeighborhoodBase<PointServerType,AccType> cn(ps,acc);

  // Initialize outputs and auxillary storage
  ps.basePoint(bp);
  cn.inflateAll();

  CALLGRIND_STOP_INSTRUMENTATION

  // Save the output
  MatlabIO::Save writer("CN_PSPreordered_output.mat");
  writer.saveScalar(pvalue,"pvalue");
  writer.save(cn.chisqNeighborhood(),"chisqN");
  writer.save(cn.chisqPoint(),"chisqP");
}

