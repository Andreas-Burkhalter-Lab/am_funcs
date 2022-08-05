#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerImageV.h"
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
  typedef Eigen::Matrix<PositionType,CoordV::RowsAtCompileTime,1> CoordVP;

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
  PointServerImageV<CoordP,CoordV> ps(gn,(ImDataType*) mxGetData(im),n_values);

  // Test tying functionality
  cout << "Max tying points: " << gn.maxTyingPoints() << endl;
  cout << gn.mSDistance << endl;

  // Do stuff with it
  sztmp[0] = 5;
  sztmp[1] = 20;
  ps.basePositionG(sztmp);
  //cout << "base point position " << ps.basePoint().position().transpose() << endl;
  //cout << "base point value " << ps.basePoint().value().transpose() << endl;
  CoordV vraw(n_values);
  CoordVP v(n_values);
  CoordVint vi(n_values);
  for (i = 0; i < 9; i++, ps++) {
    //vraw = ps.currentPoint();
    v = ps.currentPoint();
    //v = ps.currentPoint().value();
    //vi = v.cast<int>();
    vi = ps.currentPoint().cast<int>();
    cout << "Current index: " << ps.currentIndex() << ", current position: " << ps.currentPosition().transpose() << ", current value: " << vi.transpose() << endl;
  }
}

