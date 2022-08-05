/** \brief A collection of different ways to measure the distance between two points
 *
 * Metrics measure distances between two points in multidimensional space. They support the following interface:
 *   - distance(point1,point2): the actual distance
 *   - squaredDistance(point1,point2): the square distance. For some metrics (e.g., Euclidean) this can be calculated more efficiently than first calculating the distance and then squaring.
 *   - fastDistance(point1,point2): a fast value that preserves ordering, meaning that if distance(a,b) < distance(a,c), then fastDistance(a,b) < fastDistance(a,c).  For example, with Euclidean metrics, the fastDistance might be the square-distance, since it avoids the need to take a square root.  For a metric based on the L4 norm, fastDistance would be the fourth power of distance.
 *   - update(accumulator): a way of updating any parameters needed for the metric
 *
 * The representation of a point uses Eigen, e.g., Eigen::Matrix<float,n_dims,1>.
 */

#ifndef Metric_h
#define Metric_h

#include <cmath>
#include <limits>
#include "Eigen/Dense"
#include "AccumulatorMoments.h"

namespace Metric {

/** \brief The Euclidean metric, based on the L2 norm.
 * \details The distance between points \f${\bf x}\f$ and \f${\bf y}\f$ is defined as
 * \f[ d({\bf x},{\bf y}) = \sqrt{\sum_i (x_i - y_i)^2}, \f]
 * where \f$x_i\f$ is the \f$i\f$th coordinate of \f${\bf x}\f$.
 */
class Euclidean {
public:
  /** Default constructor (no action required)
   */
  Euclidean() {;}

  /** The distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The Euclidean distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return sqrt(fastDistance(x1,x2));
  }

  /** The square-distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square-distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return fastDistance(x1,x2);
  }

  /** A fast order-preserving distance (here, the square distance)
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square-distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar fastDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return (x1-x2).squaredNorm();
  }

  template <typename AccumulatorType>
  void update(const AccumulatorType& acc) { ; }
};


/** \brief The Mahalanobis distance.
 * \details The distance between points \f${\bf x}\f$ and \f${\bf y}\f$ is defined as
 * \f[ d({\bf x},{\bf y}) = \sqrt{({\bf x}-{\bf y})^T {\bf C}^{-1} ({\bf x}-{\bf y})}, \f]
 * where \f${\bf C}\f$ is the covariance matrix.  A convenient way to supply this is via AccumulatorMoments.
 */
template <typename Derived,typename Moments::CovarianceModel ctype = Moments::Full>
class Mahalanobis {
public:
  // Representation of a (vector) data point
  typedef Derived PointType;
  // The scalar, i.e., the data type for individual coordinates
  typedef typename Derived::Scalar Scalar;
  // The matrix type
  typedef typename Eigen::Matrix<Scalar,PointType::RowsAtCompileTime,PointType::RowsAtCompileTime> MatrixType;
  // The type of the covariance matrix decomposition
  typedef typename Eigen::LDLT<MatrixType> CovarType;

  /** Constructor
   * \param d The dimensionality
   */
  Mahalanobis(int d) { set(Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,Derived::RowsAtCompileTime>::Identity(d,d)); }

  /** The distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The Euclidean distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return sqrt(fastDistance(x1,x2));
  }

  /** The square-distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square-distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return fastDistance(x1,x2);
  }

  /** A fast order-preserving distance (here, the square distance)
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square-distance between x1 and x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar fastDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    mDx = x1-x2;
    mDx = mCov.transpositionsP() * mDx;
    mCov.matrixL().solveInPlace(mDx);
    return (mDx.array().square() * mID.array()).sum();
  }

  /** Update the covariance from an accumulator
   * \param acc An accumulator of the same covariance type
   */
  void update(const AccumulatorMoments<Derived,Moments::Full>& acc) { 
    mCov = acc.covariance();
    mID = acc.inverse_diagonal();
  }

  /** Set the covariance to a specified value
   * \param C The covariance (for Full, a matrix; for Diagonal, a column vector; for Isotropic, a scalar)
   */
  template <typename Derived2>
  void set(const Eigen::MatrixBase<Derived2>& C) {
    mCov.compute(C);
    mID = mCov.vectorD().cwiseInverse();
  }

private:
  mutable Eigen::Matrix<Scalar,Eigen::Dynamic,1> mDx;
  // LDLT decomposition of the covariance
  CovarType        mCov;
  // Storage for the inverse covariance
  PointType        mID;           // size d (inverse of diagonals)
};

template <typename Derived>
class Mahalanobis<Derived,Moments::Diagonal> {
public:
  // Representation of a (vector) data point
  typedef Derived PointType;
  // The scalar, i.e., the data type for individual coordinates
  typedef typename Derived::Scalar Scalar;
  typedef PointType CovarType;

  Mahalanobis(int d) { mICov = PointType::Ones(d,1); }

  template <typename D1,typename D2>
  typename D1::Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return sqrt(fastDistance(x1,x2));
  }

  template <typename D1,typename D2>
  typename D1::Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return fastDistance(x1,x2);
  }

  template <typename D1,typename D2>
  typename D1::Scalar fastDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    mDx = x1-x2;
    return (mDx.array().square() * mICov.array()).sum();
  }

  void update(const AccumulatorMoments<Derived,Moments::Diagonal>& acc) { 
    mICov = acc.inverse_diagonal();
  }

  template <typename Derived2>
  void set(const Eigen::MatrixBase<Derived2>& D) {
    mICov = D.cwiseInverse();
  }

private:
  mutable Eigen::Matrix<Scalar,Eigen::Dynamic,1> mDx;
  // Storage for the inverse covariance
  PointType        mICov;         // size d
};

template <typename Derived>
class Mahalanobis<Derived,Moments::Isotropic> {
public:
  // Representation of a (vector) data point
  typedef Derived PointType;
  // The scalar, i.e., the data type for individual coordinates
  typedef typename Derived::Scalar Scalar;
  // The type of the covariance matrix decomposition
  typedef Scalar CovarType;

  Mahalanobis() { mICov = 1.0; }
  Mahalanobis(int d) { mICov = 1.0; }

  template <typename D1,typename D2>
  typename D1::Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return sqrt(fastDistance(x1,x2));
  }

  template <typename D1,typename D2>
  typename D1::Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return fastDistance(x1,x2);
  }

  template <typename D1,typename D2>
  typename D1::Scalar fastDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return (x1-x2).squaredNorm() * mICov;
  }

  void update(const AccumulatorMoments<Derived,Moments::Isotropic>& acc) { 
    mICov = acc.inverse_diagonal();
  }

  void set(Scalar C) {
    mICov = 1.0/C;
  }

private:
  // Storage for the inverse covariance
  Scalar    mICov;
};


/** \brief A metric based on the Lp norm.
 * \details The distance \f$d\f$ between points \f${\bf x}\f$ and \f${\bf y}\f$ is defined as
 * \f[ d({\bf x},{\bf y}) = \left(\sum_i |x_i - y_i|^p\right)^{1/p}, \f]
 * where \f$x_i\f$ is the \f$i\f$th coordinate of \f${\bf x}\f$.
 * \note For \c p = 1, this is often called the "Manhattan distance."
 * \note \c p = 2 is the Euclidean distance
 * \note For \c p > 2, greater emphasis is placed on coordinates with the largest differences.  The logical extreme is \c p = \c Eigen::Infinity, which corresponds to the maximum absolute value difference.
 */
template <int p>
class LpDistance {
public:

  /** Default constructor (no action required)
   */
  LpDistance() {;}

  /** The distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The Lp norm of x1-x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return (x1-x2).template lpNorm<p>();
  }

  /** The square-distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square of the Lp norm of x1-x2
   */
  template <typename D1,typename D2>
  typename D1::Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return Eigen::internal::lpNorm_selector<D1,p>::fast2squaredNorm(fastDistance(x1,x2));
  }

  /** A fast order-preserving "distance"
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The "fast distance" between x1 and x2, typically \f$d^p\f$, except for \c p = \c Eigen::Infinity in which case it is simply \f$d\f$
   */
  template <typename D1,typename D2>
  typename D1::Scalar fastDistance(const D1 &x1,const D2 &x2) const {
    return (x1-x2).template lpNormFast<p>();
  }

  template <typename AccumulatorType>
  void update(const AccumulatorType& acc) { ; }
};



/** \class LengthScaledDirect
 * \brief A metric in which the length of each vector may be scaled
 * \details The distance between points \f${\bf x}\f$ and \f${\bf y}\f$ is defined as
 * \f[ d({\bf x},{\bf y})^p = |\hat {\bf x} - \hat{\bf y}|^p + \left[c \log(|{\bf x}|/|{\bf y}|)\right]^p \f]
 * where \f$\hat {\bf x}\f$ is the unit vector pointing in the direction of \f${\bf x}\f$, all lengths are measured using the Lp-norm, and \f$c\f$ is a scalar constant which determines the penalty for mismatched length.
 *
 * \note The distance between the "zero vector" and any other point is undefined.
 * \note For \c p = \c Eigen::Infinity, the \f$p\f$th powers are skipped and \c abs() is used instead.
 *
 * \tparam Scalar The scalar type used to represent coordinates
 * \tparam p An integer specifying which of the Lp norms to use to measure length
 */
template <typename Scalar,int p>
class LengthScaledDirect {
public:
  /** Constructor
   * \param c Specifies the size of the penalty for mismatched length.
   */
  LengthScaledDirect(Scalar c) : mC(c) {;}

  /** The distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The distance between x1 and x2
   */
  template <typename D1,typename D2>
  Scalar distance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return Eigen::internal::lpNorm_selector<D1,p>::fast2norm(fastDistance(x1,x2));
  }

  /** The square-distance between two points
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The square-distance between x1 and x2
   */
  template <typename D1,typename D2>
  Scalar squaredDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    return Eigen::internal::lpNorm_selector<D1,p>::fast2squaredNorm(fastDistance(x1,x2));
  }

  /** A fast order-preserving distance. See more complete description in LpDistance.
   * \param x1 A vector representing the first point
   * \param x2 A vector representing the second point
   * \return The "fast distance" between x1 and x2
   */
  template <typename D1,typename D2>
  Scalar fastDistance(const Eigen::MatrixBase<D1> &x1,const Eigen::MatrixBase<D2> &x2) const {
    Scalar x1norm,x2norm;

    x1norm = x1.template lpNorm<p>();
    x2norm = x2.template lpNorm<p>();
    return (x1/x1norm - x2/x2norm).template lpNormFast<p>() +
      Eigen::internal::lpNorm_selector<D1,p>::norm2fast(mC*log(x1norm/x2norm));
  }

  template <typename AccumulatorType>
  void update(const AccumulatorType& acc) { ; }

  /** Set the coefficient penalizing length mismatch
   * \param c The coefficient of the log term
   */
  void coefficient(Scalar c) {mC = c;}

private:
  Scalar mC;
};

} // namespace Metric

#endif // Metric_h

// Length-scaled metrics allow the overall length of the points to
// be adjusted to minimize their difference.  Such metrics can be
// derived from the 1-d family
//    dist(x,y)^p = Lpnorm(exp(alpha)*x-exp(-alpha)*y)^p/(Lpnorm(x)*Lpnorm(y))^(p/2) + (c*alpha)^p
// where the first term is the Lp norm of the scaled-difference, and
// the second term is the penalty applied to the scaling (c=0 means no
// penalty).
//
// In the Euclidean case, one can analytically minimize the first term
// with respect to alpha, and plug the solution for alpha into the
// penalty term.  This is the LengthScaledEuclideanDirect metric.
//
// Alternatively, one could choose alpha to minimize the total
// distance, but this requires a numerical solution for alpha and will
// therefore be slower.  One circumstance in which this is unavoidable
// is if you expect to need to have a finite distance to the zero
// vector.
//template <typename dataType,int p = 2>
//class LengthScaledLp

// Also add coordinate-scaling metrics?
