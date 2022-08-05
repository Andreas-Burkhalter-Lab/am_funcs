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

template <typename PointServer,typename Accumulator>
int clusterAllFlow(vector<int>& mapsTo,vector<int>& n,vector<typename Accumulator::Scalar>& valueCoef,PointServer& ps,Accumulator& acc,double pvalue)
{
  int n_points = ps.N();
  int nNbrs;
  int nMin = ceil(-log(pvalue));
  NeighborhoodHistory history;
  int i,j;
  CriticalNeighborhoodAlgorithms::FlowResults result;
  vector<bool> isAtMax(n_points);
  int nShells = 3;
  CriticalNeighborhoodBase<PointServer,Accumulator> cn(ps,acc);

  typename Accumulator::Scalar defaultValueCoefficient = 1;
  fill(valueCoef.begin(),valueCoef.end(),0);

  /* Test a particularly troublesome point:
  ps.basePointIndex(24784);
  ps.balanceValueCoefficient(nShells);
  nNbrs = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue,nMin);
  */

  // Move each probe point until it maps to another, or until it reaches a peak
  cout << "First pass...";
  for (i = 0; i < n_points; i++) {  // loop over probe points
    if (i%1000 == 0)
      cout << i << "..." << flush;
    ps.basePointIndex(i);
    ps.balanceValueCoefficient(nShells);
    result = CriticalNeighborhoodAlgorithms::flowToNeighbor(cn,history,pvalue,nMin);
    mapsTo[i] = cn.sortOrder()[0];
    n[i] = result.nNbrs;
    isAtMax[i] = result.isAtMax;
    if (result.isAtMax)
      valueCoef[i] = ps.valueCoefficient();
  }
  cout << "done." << endl;

  // Flow the map
  vector<int> mapsTo0 = mapsTo;   // preserve the original
  CriticalNeighborhoodAlgorithms::flowMap(mapsTo,n);
  vector<int> mapsToOld;

  // Reflow until we get consistency
  int nFlowed = 1;
  while (nFlowed > 0) {
    mapsToOld = mapsTo;
    nFlowed = 0;
    // Check points that appear to map to themselves, by flowing them all the
    // way to their peak
    for (i = 0; i < n_points; i++) {
      if (mapsTo[i] == i && !isAtMax[i]) {
	ps.basePointIndex(i);
	ps.balanceValueCoefficient(nShells);
	nNbrs = CriticalNeighborhoodAlgorithms::flowToPeak(cn,history,pvalue,nMin);
	mapsTo0[i] = cn.sortOrder()[0];
	n[i] = nNbrs;
	isAtMax[i] = true;
	valueCoef[i] = ps.valueCoefficient();
	nFlowed++;
      }
    }
    cout << nFlowed << " re-flowed to peak" << endl;
    mapsTo = mapsTo0;
    CriticalNeighborhoodAlgorithms::flowMap(mapsTo,n);
    for (i = 0; i < n_points; i++)
      if (mapsTo[i] != mapsToOld[i])
	break;
    if (i == n_points)
      break;   // map did not change even though we're still flowing points
  }

  // Return the number of points flowed all the way to the peak
  int nAtMax = 0;
  for (vector<bool>::iterator it = isAtMax.begin(); it < isAtMax.end(); it++)
    nAtMax += *it;
  return nAtMax;
}

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
  typedef PointServerImagePV<CoordP,CoordV,MetricType,true> PointServerType;
  typedef PointServerType::CoordVP CoordVP;
  typedef AccumulatorImagePV<CoordP,CoordV> AccumulatorType;

  // The p-value
  PositionType pvalue = 0.001;

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

  sz[0] = pSz[1];
  sz[1] = pSz[2];
  int n_max = sz.maxCoeff();
  int n_points = sz[0]*sz[1];
  
  // Create the PointServer
  GridNeighbors<CoordP> gn(sz,n_max);
  PointServerType ps(gn,(ImDataType*) mxGetData(im),n_values,metric);
  ps.valueCoefficient(10);  // for Euclidean
  //ps.valueCoefficient(100);  // for LengthScaled<4>
  //ps.valueCoefficient(300);  // for LengthScaled<2>
  AccumulatorType acc(2,n_values);

  vector<int> mapsTo(n_points);
  vector<int> n(n_points);
  vector<AccumulatorType::Scalar> valueCoef(n_points);

  int nAtMax = clusterAllFlow(mapsTo,n,valueCoef,ps,acc,pvalue);
  cout << nAtMax << " were flowed to their maximum" << endl;


  MatlabIO::Save writer("cluster_image_output.mat");
  writer.save(n,"n");
  writer.save(mapsTo,"map");
  writer.save(valueCoef,"valueCoef");
}

