#include <iostream>
#include "PointServerPreorderedC.h"

using namespace std;

int main()
{
  typedef float dataType;
  const int n_dims = 4;
  const int n_points = 100;
  typedef Eigen::Matrix<dataType,n_dims,Eigen::Dynamic> DataMatrixType;
  typedef PointServerPreorderedC<DataMatrixType,Eigen::VectorXi> PointServerType;

  DataMatrixType X(n_dims,n_points);
  Eigen::VectorXi sortOrder(n_points);

  //X = DataMatrixType::Random(n_dims,n_points);
  X.setRandom();
  for (int i = 0; i < n_points; i++)
    sortOrder[i] = i;

  PointServerType ps(X,sortOrder);

  cout << "The base point is " << ps.basePoint().transpose() << endl;
  PointServerType::PointType newbase = PointServerType::PointType::Zero(n_dims,1);
  ps.basePoint(newbase);
  cout << "The new base point is " << ps.basePoint().transpose() << endl;

  int n_show = 5;
  cout << "The first " << n_show << " points:\n";
  for (int i = 0; i < n_show; i++) {
    cout << ps.currentPoint().transpose() << endl;
    ps++;
  }
}

