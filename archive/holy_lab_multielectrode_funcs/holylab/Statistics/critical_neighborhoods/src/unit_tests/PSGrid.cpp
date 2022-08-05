#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "PointServerGrid.h"

//#include <valgrind/callgrind.h>

// Compile with make

using namespace std;
using namespace Eigen;

int main()
{
  // Specify the types we'll be using
  typedef float DataType;
  const int n_dims = 2;
  typedef Matrix<DataType,Dynamic,1> CoordP;
  typedef Matrix<DataType,Dynamic,Dynamic> TformType;

  int sz[] = {10,13};
  Map<Vector2i> szG(sz);
  TformType P(n_dims,n_dims);
  P <<
    1, 0,
    0, 1.8;
  GridNeighbors<CoordP> gn(szG,4,P);
  PointServerGrid<CoordP> ps(gn);

  Vector2i bp;
  bp << 4, 1;
  ps.basePointG(bp);

  for ( ; !ps.isAtEnd(); ps++) {
    cout << ps.currentIndex() << ' ' << ps.currentSquareDist() << ' ' << ps.isTiesPrevious() << endl;
    cout << ps.currentPoint().transpose() << endl;
  }
  
}

