#include <iostream>
#include "AccumulatorMoments.h"

// compile with g++ -g -I.. AccMoments.cpp -o AccMoments

int main()
{
  typedef double DataType;
  const int n_dims = 4;
  // Dynamic dimensionality
  typedef Eigen::Matrix<DataType,Eigen::Dynamic,1> PointType;
  // Static dimensionality
  //typedef Eigen::Matrix<DataType,n_dims,1> PointType;
  typedef AccumulatorMoments<PointType,Moments::Full> AccType;

  AccType acc(n_dims);
  PointType mean(n_dims),tmp(n_dims), tmp2(n_dims), tmp3(n_dims);
  DataType chisqNeighborhood,chisqPoint,dof;

  acc.restart(PointType::Zero(n_dims,1));

  mean = acc.mean();
  std::cout << "Prior to restart: " << mean.transpose() << std::endl;

  acc.restart(PointType::Constant(n_dims,2));
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;
  std::cout << "BasePoint: " << acc.basePoint().transpose() << std::endl;

  tmp = PointType::Random(n_dims);
  acc.addPoint(tmp);
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;
  tmp = PointType::Random(n_dims);
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  acc.addPoint(tmp);
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;
  tmp = PointType::Random(n_dims);
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  acc.addPoint(tmp);
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;
  tmp = PointType::Random(n_dims);
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  acc.addPoint(tmp);
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;
  tmp = PointType::Random(n_dims);
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  acc.addPoint(tmp);
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;

  tmp = PointType::Random(n_dims);
  acc.statisticsPoint(chisqPoint,dof,tmp);
  std::cout << "Before adding it, chisqPoint = " << chisqPoint << std::endl; 
  std::cout << "Adding point " << tmp.transpose() << std::endl;
  acc.addPoint(tmp);
  std::cout << acc.n() << std::endl;
  mean = acc.mean();
  std::cout << "Mean: " << mean.transpose() << std::endl;

  acc.statisticsNeighborhood(chisqNeighborhood,dof);
  std::cout << "chisqNeighborhood: " << chisqNeighborhood << ", dof: " << dof << std::endl;

  // Test dot products, etc
#if INVERT_TRUE_COVAR
  tmp2.setRandom();
  tmp3.setRandom();
  std::cout << "Random vectors:\n" << tmp2.transpose() << "\n" << tmp3.transpose() << std::endl;
  std::cout << "Their 'dot product' (relative to the data covariance matrix) = " << acc.n()*acc.dotProduct(tmp2,tmp3) << std:: endl;
  acc.solveInPlace(tmp2);
  tmp2 *= acc.n();
  std::cout << "solveInPlace (relative to data covariance matrix): " << tmp2.transpose() << std::endl;
#endif

  std::cout << "Removing point " << tmp.transpose() << std::endl;
  acc.removePoint(tmp);
  mean = acc.mean();
  std::cout << acc.n() << std::endl;
  std::cout << "Mean after removal: " << mean.transpose() << std::endl;
  acc.clear();
  mean = acc.mean();
  std::cout << "Mean after clearing: " << mean.transpose() << std::endl;
}

