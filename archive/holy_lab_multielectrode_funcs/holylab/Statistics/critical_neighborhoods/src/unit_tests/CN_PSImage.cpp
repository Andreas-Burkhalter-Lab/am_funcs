#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "AccumulatorImagePV.h"
#include "CriticalNeighborhoodBase.h"
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
  typedef PointServerImagePV<CoordP,CoordV,Metric::Euclidean> PointServerType;
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
  Metric::Euclidean metric;

  sz[0] = pSz[1];
  sz[1] = pSz[2];
  
  // Create the PointServer
  GridNeighbors<CoordP> gn(sz,10);
  PointServerType ps(gn,(ImDataType*) mxGetData(im),n_values,metric);
  ps.valueCoefficient(1);
  AccumulatorType acc(2,n_values);
  CriticalNeighborhoodBase<PointServerType,AccumulatorType> cn(ps,acc,false);

  // Now set the base point at different locations and find the 20 closest neighbors
  int n = sz[0]*sz[1];
  //int n = 100;
  std::vector<int> n_nbrs(n);
  std::vector<int> n_checked(n);
  typedef Eigen::Matrix<PositionType,2,Eigen::Dynamic> MeanPType;
  typedef Eigen::Matrix<PositionType,Eigen::Dynamic,Eigen::Dynamic> MeanVType;
  MeanPType meanP(2,n);
  MeanVType meanV(n_values,n);
  ImagePoint<MeanPType,MeanVType> mean(meanP,meanV);
  cout << n << endl;
  clock_t t0,t1;
  t0 = clock();
  for (i = 0; i < n; i++) {
    ps.basePointIndex(i);
    n_nbrs[i] = cn.inflateToCriterion(pvalue);
    n_checked[i] = cn.sortOrder().size();
    //mean = acc.mean();
    //meanP.col(i) = mean.position();
    //meanV.col(i) = mean.value();
    meanP.col(i) = acc.mean().position();
    meanV.col(i) = acc.mean().value();
  }
  t1 = clock();
  cout << t0 << ' ' << t1 << ' ' << CLOCKS_PER_SEC << endl;
  cout << "Total time: " << (int)((t1-t0) / (CLOCKS_PER_SEC/1000)) << "ms" << endl;

  MatlabIO::Save writer("CN_Image_output.mat");
  writer.save(n_nbrs,"n_nbrs");
  writer.save(n_checked,"n_checked");
  writer.save(meanP,"meanP");
  writer.save(meanV,"meanV");
}

