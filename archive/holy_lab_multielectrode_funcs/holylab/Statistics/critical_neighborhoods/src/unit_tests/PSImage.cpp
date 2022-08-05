#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerImagePV.h"
#include "Metric.h"
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
  PointServerImagePV<CoordP,CoordV,Metric::Euclidean> ps(gn,(ImDataType*) mxGetData(im),n_values,metric);

  // Do stuff with it
  sztmp[0] = 5;
  sztmp[1] = 20;
  ps.basePointG(sztmp);
  ps.valueCoefficient(0.5);
  //cout << "base point position " << ps.basePoint().position().transpose() << endl;
  //cout << "base point value " << ps.basePoint().value().transpose() << endl;
  CoordV v(n_values);
  CoordVint vi(n_values);
  for (i = 0; i < 5; i++, ps++) {
    //v = ps.currentPoint().value();
    //vi = v.cast<int>();
    vi = ps.currentPoint().value().cast<int>();
    cout << "Current index: " << ps.currentIndex() << ", current position: " << ps.currentPoint().position().transpose() << ", current value: " << vi.transpose() << endl;
  }

  // Now set the base point at different locations and find the 20 closest neighbors
  int n = sz[0]*sz[1];
  //int n = 1000;
  std::vector<int> nbrIndex(n);
  cout << n << endl;
  clock_t t0,t1;
  t0 = clock();
CALLGRIND_START_INSTRUMENTATION
  for (i = 0; i < n; i++) {
    ps.basePointIndex(i);
    for (int j = 0; j < 20; j++)
      ps++;
    nbrIndex[i] = ps.currentIndex();
  }
CALLGRIND_STOP_INSTRUMENTATION
  t1 = clock();
  cout << t0 << ' ' << t1 << ' ' << CLOCKS_PER_SEC << endl;
  cout << "Total time: " << (int)((t1-t0) / (CLOCKS_PER_SEC/1000)) << "ms" << endl;

  MatlabIO::Save writer("CN_Image_output.mat");
  writer.save(nbrIndex,"nbrIndex");
}
