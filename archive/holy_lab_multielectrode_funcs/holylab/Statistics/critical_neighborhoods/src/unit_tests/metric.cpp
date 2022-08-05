#include <iostream>
#include "Metric.h"
#include "AccumulatorMoments.h"

// Compile with
//   g++ -g -I.. metric.cpp -o metric

using namespace Eigen;
using namespace std;

int main()
{
  typedef float Scalar;
  typedef Matrix<Scalar,Dynamic,1> PointType;
  const int n_dims = 3;

  PointType x1(n_dims),x2(n_dims);
  x1.setRandom();
  x2.setRandom();

  cout << "x1: " << x1.transpose() << endl;
  cout << "x2: " << x2.transpose() << endl;

  Metric::Euclidean euclidean;
  Metric::LpDistance<4> l4;
  Metric::LengthScaledDirect<Scalar,2> ls(10);
  Metric::Mahalanobis<PointType,Moments::Full> lm(n_dims);

  // Mahalanobis Full
  Matrix<Scalar,Dynamic,Dynamic> xr(n_dims,10);
  xr.setRandom();
  PointType xrmean = xr.rowwise().sum()/xr.cols();
  xr.colwise() -= xrmean;
  Matrix<Scalar,Dynamic,Dynamic> C = xr*xr.transpose();
  cout << "Covariance matrix (Mahalanobis):" << endl
       << C << endl;
  lm.set(C);

  /*
  // Mahalanobis Diagonal
  Matrix<Scalar,Dynamic,1> md(n_dims);
  md.setRandom();
  md = md.cwiseAbs();
  cout << "Diagonal of covariance matrix (Mahalanobis):" << endl
       << md.transpose() << endl;
  lm.set(md);
  */

  /*
  // Mahalanobis Isotropic
  lm.set(3);
  */

  cout << "Euclidean dist: " << euclidean.distance(x1,x2) << endl;
  cout << "Mahalanobis squareddist: " << lm.squaredDistance(x1,x2) << endl;
  cout << "L4 fast: " << l4.fastDistance(x1,x2) << endl;
  cout << "L4 distance: " << l4.distance(x1,x2) << endl;
  cout << "L4 squaredDistance: " << l4.squaredDistance(x1,x2) << endl;
  cout << "LengthScaledEuclideanDirect: " << ls.distance(x1,x2) << endl;
  cout << "Comparison: " << (x1/x1.norm() - x2/x2.norm()).norm() << ' ' << 10*log(x1.norm()/x2.norm()) << endl;
  x2.setZero();
  cout << "LengthScaledEuclideanDirect, x2 zero: " << ls.distance(x1,x2) << endl;

  /*  
  Vector2f x,y;
  Metric::Saturating<PointType> msat(euclidean,0,1,2);
  x[0] = 1;
  x[1] = 0.5;
  y[0] = 0.8;
  y[1] = 0.3;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  y[0] = 0.3;
  y[1] = 0.8;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  y[1] = 1;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  x[1] = 1;
  y[0] = 1;
  y[1] = 0.8;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  x[0] = 0;
  x[1] = 0.5;
  y[0] = 0.8;
  y[1] = 0.3;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  y[0] = 0.3;
  y[1] = 0.8;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  y[1] = 0;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  x[1] = 0;
  y[0] = 0;
  y[1] = 0.8;
  cout << "x: " << x.transpose() << ", y: " << y.transpose() << endl;
  cout << msat.distance(x,y) << endl;
  */
}
