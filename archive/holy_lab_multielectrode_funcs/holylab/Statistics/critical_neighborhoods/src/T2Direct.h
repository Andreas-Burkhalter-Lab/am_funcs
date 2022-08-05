/**  Directly calculate T^2 and moments from a set of points
 *
 * #include "T2Direct.h" <BR>
 *
 * \note It is not recommended to use these classes with integer data types due to the problems with overflow.
 *  
 */

#ifndef T2Direct_h
#define T2Direct_h

// SYSTEM INCLUDES
//
#include "Eigen/Dense"
#include "Eigen/Cholesky"

// PROJECT INCLUDES
//
#include "AccumulatorMoments.h"

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//


/// \brief Calculate the mean, covariance, and T^2 for a set of points

/// \tparam Derived The data matrix type.  This must manage its own storage; do not supply a Map type for this, even if the points you pass to addPoint will be of type Map<Derived>.  Do not use integer types as the scalar (use float/double instead).
/// \tparam ctype The enum type used to indicate the covariance model (see namespace Moments)
template <typename Derived,Moments::CovarianceModel ctype = Moments::Isotropic>
class T2Direct
{
public:
// Typedefs and enums
  /** The scalar, i.e., the data type for individual coordinates  */
  typedef typename Derived::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,1> VectorType;
  typedef Scalar CovarianceType; 

// LIFECYCLE
  /** Default constructor (not available with "Full") */
  T2Direct() {;}

  /** Constructor specifying dimensionality
   * \param d The number of coordinates in the representation of a point (the embedding dimensionality)
   * \param maxpoints The largest number of points in any data set
   */
  T2Direct(int d,int maxpoints) {;}
  
// INQUIRY
  int nMin() { return 1; }

// OPERATORS

// OPERATIONS
  // We template it so that it works with Eigen::Block types, for example
  template <typename OtherDerived>
  Scalar calculate(const OtherDerived &x, VectorType &mu, CovarianceType &C) {
    n = x.cols();
    mu = x.rowwise().sum()/n;
    C = (x.colwise() - mu).array().square().sum()/(x.size()-1);
    return n*(mu.array().square().sum()/C);
  }
  Scalar calculate_newmu(const VectorType &mu,const CovarianceType &C) {
    return n*(mu.array().square().sum()/C);
  }

private:
  int n;
};


// Make sure Doxygen skips all the template specializations
#ifndef DOXYGEN_SHOULD_SKIP_THIS

//
// Template specialization for Diagonal covariance
//
template <typename Derived>
class T2Direct<Derived,Moments::Diagonal>
{
public:
// Typedefs and enums
  typedef typename Derived::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,1> VectorType;
  typedef VectorType CovarianceType;

// LIFECYCLE
  T2Direct() {;}
  T2Direct(int d,int maxpoints) {;}
  
// INQUIRY
  int nMin() { return 1; }

// OPERATORS

// OPERATIONS
  template <typename OtherDerived>
  Scalar calculate(const OtherDerived &x, VectorType &mu, CovarianceType &C) {
    n = x.cols();
    mu = x.rowwise().sum()/n;
    C = (x.colwise() - mu).array().square().rowwise().sum()/(n-1);
    return n*(mu.array().square()/C.array()).sum();
  }
  Scalar calculate_newmu(const VectorType &mu,const CovarianceType &C) {
    return n*(mu.array().square()/C.array()).sum();
  }

private:
  int n;
};


//
// Template specialization for Full covariance
//
template <typename Derived>
class T2Direct<Derived,Moments::Full>
{
public:
// Typedefs and enums
  typedef typename Derived::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,1> VectorType;
  typedef typename Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,Derived::RowsAtCompileTime> CovarianceType;

// LIFECYCLE
  T2Direct(int d,int maxpoints) : mDx(d,maxpoints),mCinvMu(d,1),mCDecomp(d) {;}
  
// INQUIRY
  int nMin() { return mDx.size(); }

// OPERATORS

// OPERATIONS
  template <typename OtherDerived>
  Scalar calculate(const OtherDerived &x, VectorType &mu, CovarianceType &C) {
    n = x.cols();
    mu = x.rowwise().sum()/n;
    mDx.leftCols(n) = x.colwise()-mu;
    C = (mDx.leftCols(x.cols())*mDx.leftCols(n).adjoint())/(n-1);
    mCDecomp.compute(C);
    mCinvMu = mCDecomp.solve(mu);
    return x.cols()*mu.dot(mCinvMu);
  }
  Scalar calculate_newmu(const VectorType &mu, const CovarianceType &C) {
    // This doesn't use C
    mCinvMu = mCDecomp.solve(mu);
    return n*mu.dot(mCinvMu);
  }
private:
  int        n;
  Derived    mDx;
  VectorType mCinvMu;
  Eigen::LDLT<CovarianceType> mCDecomp;
};



#endif // DOXYGEN_SHOULD_SKIP_THIS  (for template specializations)


#endif  // T2Direct_h
