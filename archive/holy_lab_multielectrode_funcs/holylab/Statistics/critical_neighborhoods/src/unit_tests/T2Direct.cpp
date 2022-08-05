#include <iostream>
#include "T2Direct.h"

using namespace std;
using namespace Eigen;

int main() {
  typedef double Scalar;
  const int n_dims = 5;
  Matrix<Scalar,Dynamic,Dynamic> x(n_dims,30);
  Matrix<Scalar,Dynamic,1>       mu(n_dims,1);
  //Matrix<Scalar,Dynamic,1> C(n_dims,1);
  Scalar C;
  T2Direct<Matrix<Scalar,Dynamic,Dynamic>,Moments::Isotropic> T2D(n_dims,100);
  Scalar T2;

  x.setRandom();
  T2 = T2D.calculate(x,mu,C);
  cout << "mu = " << mu.transpose() << endl;
  cout << "C = " << C << endl;
  cout << "T2 = " << T2 << endl;

  for (int i = 0; i < 100; i++) {
    x.setRandom();
    T2 = T2D.calculate(x,mu,C);
    cout << "T2 = " << T2 << endl;
  }
}
  
