#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "AccumulatorImagePV.h"
#include "CriticalNeighborhoodBase.h"
#include "NeighborhoodHistory.h"
#include "mat.h"
#include "MatlabIO.h"

#include <valgrind/callgrind.h>
#include <time.h>

// Compile with make

using namespace std;

int main()
{
  // Specify the types we'll be using
  typedef uint8_t ImDataType;
  typedef float   PositionType;
  typedef Eigen::Matrix<PositionType,2,1> CoordP;
  typedef Eigen::Matrix<ImDataType,Eigen::Dynamic,1> CoordV;
  typedef Eigen::Matrix<int,CoordV::RowsAtCompileTime,1> CoordVint;
  //typedef Metric::LengthScaledDirect<PositionType,2> MetricType;
  typedef Metric::Euclidean MetricType;
  typedef PointServerImagePV<CoordP,CoordV,MetricType> PointServerType;
  typedef PointServerType::CoordVP CoordVP;
  typedef AccumulatorImagePV<CoordP,CoordV> AccumulatorType;

  // The p-value
  PositionType pvalue = 0.001;
  int n_min = ceil(-log(pvalue));  // stability under bootstrap
  // The Manhattan-distance threshold for setting the flow
  const int displacement_thresh = 2;   // 2-pixel displacement to terminate flow

  // Load the data
  MatlabIO::Load reader("CN_Image_input.mat");
  mxArray *im = reader.load("im");

  // Prepare all the other information
  const int *pSz = mxGetDimensions(im);
  int n_values = pSz[0];
  Eigen::Vector2i sz,sztmp;
  int i;
  //MetricType metric(0.1);
  MetricType metric;
  NeighborhoodHistory history;

  sz[0] = pSz[1];
  sz[1] = pSz[2];
  int n_max = sz.maxCoeff();
  
  // Create the PointServer
  GridNeighbors<CoordP> gn(sz,n_max);
  PointServerType ps(gn,(ImDataType*) mxGetData(im),n_values,metric);
  ps.valueCoefficient(1);  // for Euclidean
  //ps.ValueCoefficient(100);  // for LengthScaled<4>
  //ps.ValueCoefficient(300);  // for LengthScaled<2>
  AccumulatorType acc(2,n_values);
  CriticalNeighborhoodBase<PointServerType,AccumulatorType> cn(ps,acc);

  // Cycle over pixels and determine where they flow, by iterating mean shift
  //int n = sz[0]*sz[1];
  int n = 10000;
  std::vector<int> n_nbrs(n);
  int n_nbrs_tmp;
  std::vector<int> target(n);
  std::vector<int8_t> iscycling(n);
  CoordP pixelStart;

  cout << n << endl;
  bool isAtMax;
  clock_t t0,t1;
  t0 = clock();
  CALLGRIND_START_INSTRUMENTATION;
  for (i = 0; i < n; i++) {
    // Start from a new pixel
    cout << "i = " << i << ", n =";
    ps.basePointIndex(i);
    pixelStart = ps.basePointG().cast<PositionType>();  // save the grid coordinates
    //cout << "pixelStart " << pixelStart.transpose() << endl;
    history.restart();
    while (true) {
      // Run the critical neighborhood algorithm
      n_nbrs_tmp = cn.inflateToCriterion(pvalue,n_min);
      //if (n_nbrs_tmp < n_min)
      //  n_nbrs_tmp = n_min;
      cout << ' ' << n_nbrs_tmp << '(' << cn.sortOrder().size()+ps.heapSize() << ')';
      //cout << "Prior to backtracking:\n";
      //acc.print();
      // Backtrack
      cn.backtrack(n_nbrs_tmp);
      //cout << "Final state:\n";
      //acc.print();
      //cout << "Mean position: " << acc.mean().position().transpose() << endl;
      //cout << "Mean value: " << acc.mean().value().transpose() << endl;
      // Determine whether we are cycling and at the peak
      history.add(cn.sortOrder(),n_nbrs_tmp);
      isAtMax = history.isAtMax();
      if (isAtMax)
	break;
      // Set the base point to the mean of the neighborhood
      ps.basePoint(acc.mean());
      // Determine whether we have gone far enough
      if ( (ps.basePoint().position()-pixelStart).array().abs().sum() >= displacement_thresh )
	break;
    }
    // Set the flow information
    target[i] = ps.basePointIndex();
    n_nbrs[i] = n_nbrs_tmp;
    iscycling[i] = isAtMax;
    cout << endl;
  }
  CALLGRIND_STOP_INSTRUMENTATION;
  t1 = clock();
  cout << t0 << ' ' << t1 << ' ' << CLOCKS_PER_SEC << endl;
  cout << "Total time: " << (int)((t1-t0) / (CLOCKS_PER_SEC/1000)) << "ms" << endl;

  MatlabIO::Save writer("flow_image_output.mat");
  writer.save(n_nbrs,"n_nbrs");
  writer.save(target,"map0");
  writer.save(iscycling,"iscycling");

  // Now flow the map, and determine which point to themselves
}

