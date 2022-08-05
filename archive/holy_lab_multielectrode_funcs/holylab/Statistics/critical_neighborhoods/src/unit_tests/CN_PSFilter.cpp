#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerFilter1d.h"
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
  typedef double DataType;
  typedef MatlabIO::MapMx<DataType>::MapMatrixType DataMatrixType;
  typedef MatlabIO::MapMx<DataType>::ColVectorType PointType;
  typedef Metric::Euclidean                        MetricType;
  typedef PointServerFilter1d<DataMatrixType> PointServerType;
  const Moments::CovarianceModel cmodel = Moments::Isotropic;
  typedef AccumulatorMoments<PointType,cmodel> AccType;

  // Specify the p-value
  const DataType pvalue = 0.001;

  // Load the data
  DataMatrixType X(NULL,0,0);
  MatlabIO::Load reader("PSFilter_input.mat");
  reader.load("x",X);
  int n_dims = X.rows();
  int n_points = X.cols();

  // Initialize the variables needed for critical neighborhood analysis
  PointServerType ps(X);
  AccType acc(n_dims);
  NeighborhoodStatistics<DataType> ns(pvalue,cmodel);
  CriticalNeighborhoodBase<PointServerType,AccType> cn(ps,acc,ns);
  NeighborhoodHistory history;

  // Initialize outputs and auxillary storage
  Eigen::Matrix<DataType,1,Eigen::Dynamic> Xfiltered(1,n_points);
  Eigen::VectorXi n_checked(n_points);
  Eigen::VectorXi n(n_points);
  Eigen::VectorXi n_iter(n_points);

  //CALLGRIND_START_INSTRUMENTATION

  // Satisfy the critical neighborhood criterion for each point
  for (int i = 0; i < n_points; i++) {
    //cout << "i = " << i << ":" << endl;
    ps.basePointIndex(i);
    n_checked[i] = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history);
    n[i] = cn.backtrackPeak();
    Xfiltered[i] = acc.mean()[0];
    n_iter[i] = history.n_history().size();
  }

  //CALLGRIND_STOP_INSTRUMENTATION;


  // Save the output
  MatlabIO::Save writer("CN_PSFilter_output.mat");
  writer.save(n_checked,"n_checked");
  writer.save(n,"n");
  writer.save(n_iter,"n_iter");
  writer.saveScalar(pvalue,"pvalue");
  writer.save(Xfiltered,"xfiltered");
}

