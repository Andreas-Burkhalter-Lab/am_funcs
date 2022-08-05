/**  Accumulate the mean and covariance of the displacement
 *
 * #include "AccumulatorMoments.h" <BR>
 *
 * \note It's not recommended to use AccumulatorMoments with integer
 * data types, because of overflow.  Define the accumulator as type
 * float or double, and (using a wrapper class) cast the individual
 * data points as they are added.
 *  
 */

#ifndef AccumulatorMoments_h
#define AccumulatorMoments_h

// The following flag determines whether pointwise statistics use the
// "sample" covariance matrix (mean-centered) or the
// basepoint-centered covariance. Basepoint-centered covariance is
// preferred, because it's both more robust and exploits the notion
// that "true center" is best determined from many points (e.g., the
// previous iteration of mean shift).
// The main negative is that, for diagonal covariance models, the
// inverse covariance needs to be computed twice.  But since it is
// diagonal, this is not terribly expensive.
#define PTWISESTAT_MEANCENTERED 0

// SYSTEM INCLUDES
//
#include <limits>
#include "Eigen/Dense"

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//

/// \namespace Moments Types used to accumulate a model of a neighborhood based on moments
namespace Moments {

  /// Flags that indicate the nature of the covariance model to be used.  In reality all of the covariances listed below are scaled by \f$ 1 + {1\over n^2}\f$, arising from the asymptotic expansion of \f$\left(1-{1\over n}\right) k(n)\f$.
enum CovarianceModel {
  Isotropic, ///< Forces the covariance matrix to be proportional to the identity.  In this case, \f[T^2 = {1\over \sigma^2} \left(\sum_i \Delta {\bf x}_i\right)^2,\f] with \f[\sigma^2 = {1\over d} \sum_i \Delta {\bf x}_i^2.\f]   \f$\Delta {\bf x}_i\f$ is the vector displacement of the \f$i\f$th point from the base point of the neighborhood, and the square of a vector is its dot-product with itself.
  Diagonal,  ///< Off-diagonal elements of the covariance matrix are forced to be zero, but allows differences among the diagonal elements.  If \f$\Delta x_{ij}\f$ is the \f$j\f$th coordinate of the \f$i\f$th point's displacement, then \f[T^2 = \sum_j {\left(\sum_i \Delta x_{ij}\right)^2 \over \sigma_j^2},\f] where\f[\sigma_j^2 = \sum_i \Delta x_{ij}^2.\f]  This makes the total \f$T^2\f$ equal to the sum of the \f$\chi_j^2\f$ for each coordinate.\n For AccumulatorMoments, one can optionally specify that certain coordinates should be "grouped," modeling the covariances as isotropic among subsets of coordinates.  For example, if the first two coordinates are grouped, and all of the rest each define their own group (i.e., they are not grouped at all), then \f[T^2 = {\left(\sum_i \Delta x_{i1}\right)^2 + \left(\sum_i \Delta x_{i2}\right)^2 \over {1\over 2}\sum_i \left(\Delta x_{i1}^2 + \Delta x_{i2}^2\right)} +  \sum_{j>2} {\left(\sum_i \Delta x_{ij}\right)^2 \over \sum_i \Delta x_{ij}^2}.\f]  If all coordinates are placed in a single group, this is equivalent to the Isotropic case.\n Grouping can be useful if, for example, your first two coordinates correspond to the \f$x,y\f$ position in space and the third corresponds to some measurement (e.g., image brightness) made at that position---this gives you the flexibility to say that \f$x-\f$ and \f$y-\f$displacements are to be treated equivalently, but that spatial distance is in a different category from brightness.  Typically, groups will correspond to sets of coordinates measured in common units.
    Full,      ///< Allow a full covariance matrix, \f[T^2 = \Delta{\bf X}^T {\bf S}^{-1} \Delta{\bf X},\f] where \f[\Delta{\bf X} = \sum_{i=1}^n \Delta{\bf x}_i\f] and \f[{\bf S} = \sum_{i=1}^n \Delta{\bf x}_i\otimes\Delta{\bf x}_i.\f] This is the most complete model, but it can be significantly slower and requires the most data.  It's also worth knowing that if \f$n<d\f$ (\f$d\f$ being the dimensionality) vectors \f$\Delta {\bf x}_i\f$ have been accumulated, and we use the Moore-Penrose pseudoinverse for calculating \f${\bf S}^{-1}\f$, then \f$T^2 = n\f$ (before the bootstrap correction).  Consequently, in general the critical neighborhood criterion cannot be triggered until \f$n>d\f$.  For this and other reasons (memory requirements and speed), the Full covariance model is recommended only for spaces of modest dimensionality; for high-dimensional spaces you may prefer Diagonal or FullRegularized.
  FullRegularized ///< ? not yet implemented. Consider including sparseness, so that it might work for high-dimensional spaces? Regularize either by mixing in a component of the diagonal or by limiting to a few singular values.
};

}  // namespace Moments


/// \brief Accumulate the mean and covariance of a set of points

/// \tparam Derived The column vector type, e.g., Vector4f.  This must manage its own storage; do not supply a Map type for this, even if the points you pass to addPoint will be of type Map<Derived>.  Do not use integer types as the scalar (use float/double instead).
/// \tparam ctype The enum type used to indicate the covariance model (see namespace Moments)
template <typename Derived,Moments::CovarianceModel ctype = Moments::Isotropic>
class AccumulatorMoments
{
public:
// Typedefs and enums
  /** Representation of a (vector) data point */
  typedef Derived PointType;

  /** The scalar, i.e., the data type for individual coordinates  */
  typedef typename Derived::Scalar Scalar;

  /** A dynamically-sized vector, even if the class is defined using fixed size types */
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> DynamicPointType;  // not needed for Isotropic

  /** A column vector of integers, used to specify groups */
  typedef Eigen::Matrix<int,Derived::RowsAtCompileTime,1> ColumnOfInts;

  /** The return type of the expression for the mean */
  typedef typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<Scalar>, const PointType, const Eigen::CwiseUnaryOp<Eigen::internal::scalar_quotient1_op<Scalar>, const PointType> > MeanXpr;

  /** The return type of the second moment */
  typedef Scalar CovarianceType;

// LIFECYCLE

  /** Default constructor
   *  This is the "typical" constructor if you know the dimensionality at compile time.
   *  \param lambdaRatioThresh Any eigenvalues/eigenvectors of the covariance matrix satisfying \f$\lambda_i/\lambda_\mathrm{max} < \f$ lambdaRatioThresh will be excluded (default value 0, i.e., use all non-zero eigenvalues/eigenvectors).  Choosing a value bigger than zero may only make sense if your observations all have the same scale (and/or same units) so that the comparison is meaningful.  For Isotropic, lambdaRatioThresh is not used, but (for consistency) the constructor still accepts it as a valid parameter.
   */
  AccumulatorMoments(Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) { clear(); }  // for Isotropic lambdaRatioThresh is unused

  /** Constructor specifying dimensionality
   *  This is the "typical" constructor when the dimensionality is specified at run time.
   * \param d The number of coordinates in the representation of a point (the embedding dimensionality)
   *  \param lambdaRatioThresh Any eigenvalues/eigenvectors of the covariance matrix satisfying \f$\lambda_i/\lambda_\mathrm{max} < \f$ lambdaRatioThresh will be excluded (default value 0).
   */
  AccumulatorMoments(int d,Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) : mBasePoint(d), mDx(d), mDxCum(d) { clear(); }
  
  /** A more "advanced" constructor specifying coordinate grouping
   *  \param group A vector of length d (= dimensionality), containing integers which index isotropic groups.
   *  \param lambdaRatioThresh Any eigenvalues/eigenvectors of the covariance matrix satisfying \f$\lambda_i/\lambda_\mathrm{max} < \f$ lambdaRatioThresh will be excluded (default value 0).
   *  \details For example,
   *  \code
   *    Eigen::VectorXi group(5);
   *    group << 0, 1, 1, 2, 1;
   *    AccumulatorMoments<float,Eigen::Dynamic,Moments::Diagonal> acc(group);
   *  \endcode
   *  would model the covariance matrix as diag([a b b c b]) (Matlab notation), where b is the mean square displacement of the 2nd, 3rd, and 5th coordinates.\n To minimize storage requirements, group labels should start with 0 and not skip any values.
   *  \note Grouping is available only for CovarianceModel=Diagonal.  Note that CovarianceModel=Isotropic is equivalent to using a single group for all coordinates.
   *  \note This does not influence dof(); in the example above, for dofFlag=Embedding, dof() would still report 5.
   *  \note Presumably one cannot mix a fixed-size group with a dynamic-size acc, or vice versa.
   */
  AccumulatorMoments(const ColumnOfInts &group,Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) : mBasePoint(group.size()), mDx(group.size()), mDxCum(group.size()) { clear(); }

// OPERATORS

// OPERATIONS                       
  /** Restart from the current basepoint */
  void clear() {
    mN = 0;
    mDxCum = PointType::Zero(mDxCum.size());
    mDx2Cum = 0;
    mIsValid = false;
  }

  /** Restart with a new basepoint
   *
   * \param coords The coordinates of the new base point
   */
  template <typename T> void restart(const Eigen::MatrixBase<T>& coords) {
    clear();
    mBasePoint = coords;
  }

  /** Add a new point to the neighborhood, and incorporate into moments
   *
   * \param p A reference to the point
   */
  template <typename T> void addPoint(const Eigen::MatrixBase<T>& p) {
    mDx = p-mBasePoint;
    mN++;
    mDxCum += mDx;
    mDx2Cum += mDx.squaredNorm();
    mIsValid = false;
    //std::cout << "n = " << mN << ", mDx " << mDx.transpose() << ", mDxCum " << mDxCum.transpose() << ", mDx2Cum " << mDx2Cum << std::endl;
  }

  /** Add a new point, allowing for a scalar weight
   *
   * \param p A reference to the point
   * \param w The weight applied to the point. +1 would be a "normal"
   * point, but you could also pass the multiplicity m for a set of
   * points sharing the same location.  -1 is a point that repels,
   * rather than attracts, the center of mass.
   * \note Adding a point with weight -1 is the same as calling
   * removePoint.
   * \sa addPoint
   */
  template <typename T> void addPoint(const Eigen::MatrixBase<T>& p,Scalar w) {
    mDx = p-mBasePoint;
    mN += w;
    mDxCum += w*mDx;
    mDx2Cum += w*mDx.squaredNorm();
    mIsValid = false;
  }


  /** Remove a point from the neighborhood, subtracting its effect on moments
   *
   * \param p A reference to the point
   */
  template <typename T> void removePoint(const Eigen::MatrixBase<T>& p) {
    mDx = p-mBasePoint;
    mN--;
    mDxCum -= mDx;
    mDx2Cum -= mDx.squaredNorm();
    mIsValid = false;
    //std::cout << "n = " << mN << ", mDx " << mDx.transpose() << ", mDxCum " << mDxCum.transpose() << ", mDx2Cum " << mDx2Cum << std::endl;
  }

  /** Remove a point, allowing for a scalar weight
   * \param p A reference to the point
   * \param w The weight applied to the point
   * \sa addPoint
   */
  template <typename T> void removePoint(const Eigen::MatrixBase<T>& p,Scalar w) {
    mDx = p-mBasePoint;
    mN -= w;
    mDxCum -= w*mDx;
    mDx2Cum -= w*mDx.squaredNorm();
    mIsValid = false;
  }


  /** The covariance-weighted dot product for different left and right vectors
   * \param l the "left" vector (order does not matter, this is symmetric)
   * \param r the "right" vector
   * \return the output dot product, \f$p = {\bf l}^t {\bf C_\mathrm{cum}}^{-1} {\bf r}\f$.
   * \note \f$C_\mathrm{cum}\f$ is the "cumulative covariance" of the data points; to get the dot product with respect to the standard data covariance, multiply the result by n().
   */
  template <typename Tl,typename Tr>
  Scalar dotProduct(const Eigen::MatrixBase<Tl>& l,const Eigen::MatrixBase<Tr>& r) {
    if (!mIsValid)
      prepareInverse();
    return l.dot(r)*mICov;
  }

  /** Precalculate \f${\bf C_\mathrm{cum}}^{-1} {\bf v}\f$.
   * This is preferable if you are taking many dot products against the same vector---after performing this calculation, you can just use dot().
   * \param v the vector (note the value is modified)
   * \note \f$C_\mathrm{cum}\f$ is the "cumulative covariance" of the data points; to get the dot product with respect to the standard data covariance, multiply the result by n().
   */
  template <typename T>
  void solveInPlace(Eigen::MatrixBase<T>& v) {
    if (!mIsValid)
      prepareInverse();
    v *= mICov;
  }

  /** The covariance-weighted dot product when the "left" and "right" vectors are the same
   * \param v the vector
   * \return the output dot product, \f$p = {\bf v}^t {\bf C_\mathrm{cum}}^{-1} {\bf v}\f$.
   * \note \f$C_\mathrm{cum}\f$ is the "cumulative covariance" of the data points; to get the dot product with respect to the standard data covariance, multiply the result by n().
   */
  template <typename T>
  Scalar dotProduct(const Eigen::MatrixBase<T>& v) {
    if (!mIsValid)
      prepareInverse();
    return v.squaredNorm()*mICov;
  }


// ACCESS
  /** \return The number of points accumulated so far */
  Scalar n() const { return mN; }

  /** \return The base point */
  const PointType& basePoint() const { return mBasePoint; }

  /** \return The number of degrees of freedom */
  int dof() { if (!mIsValid) prepareInverse(); return mDof; }

  Scalar inverse_diagonal() const { return mICov; }

#if 0
  /** \return The cumulative "covariance" (the covariance times n()).  This is basePoint-centered, not mean-centered. */
  const CovarianceType covariance_n() const { return mDx2Cum; }
#endif

// INQUIRY
  /** Calculate the statistical signficance of the centroid displacement
   *
   * \param T2 This \f$T^2\f$ tests whether the center of mass of all points in the neighborhood can be distinguished from the base point
   * \param dof The number of degrees of freedom
   */
  void statisticsNeighborhood(Scalar& T2,Scalar &dof) {
    if (!mIsValid)
      prepareInverse();
    T2 = mICov * mDxCum.squaredNorm();
    dof = mDof;
  }

  /** Tests whether a point is an outlier relative to the accumulated points
   * \param T2 This \f$T^2\f$ tests whether the data point p is an outlier by the standards of the data covariance matrix. Note that while the central limit theorem guarantees that the "neighborhood statistic" (for the centroid) asymptotically obeys gaussian statistics, there is no such guarantee for pointwise statistics. So take this output with a grain of salt.
   * \param dof The number of degrees of freedom
   * \param p The point to be tested
   */
  template <typename T>
  void statisticsPoint(Scalar& T2,Scalar &dof,const Eigen::MatrixBase<T> &p) {
    if (mN == 0) {
      T2 = 0;
      dof = 0;
      return;
    }
#if PTWISESTAT_MEANCENTERED
    if (!mIsValid)
      prepareInverse();
    mDx = p - mean();
    T2 = mN * mICov * mDx.squaredNorm();  // factor of N: "std" rather than "sem"
#else
    mDx = p - mBasePoint;
    T2 = mDof * mN * mDx.squaredNorm() / mDx2Cum;
#endif
    dof = mDof;
  }

  /** Calculate the center of mass of all points accumulated thus far
   *
   * \return An expression for the mean
   */
  MeanXpr mean() const {
    // Note: when mN == 0, this returns NaN. This is actually desirable,
    // because the basepoint is not a member of the neighborhood, so in
    // fact one posesses no information about the mean.
    // The alternative would be to check for mN == 0, and in that case
    // just return the basepoint.  But this will induce very weird
    // behavior if, for example, the basepoint is off somewhere in
    // lala-land; when there are no points, it will return the
    // basepoint, but as soon as you add a single point it will "jump"
    // straight to that new point.  This would be very confusing.
    //std::cout << "mean(): basepoint " << mBasePoint.transpose() << std::endl;
    //std::cout << "mean(): mDxCum " << mDxCum.transpose() << std::endl;
    //std::cout << "mean(): mean " << (mBasePoint+mDxCum/mN).transpose() << std::endl;
    return mBasePoint+(mDxCum/mN);
  }

  /** Calculate a "scalar variance," the mean square displacement across coordinates
   *
   * \return The mean square displacement
   */
  Scalar meanSquareDisplacement() const { return mDx2Cum/mN; }
    
protected:
private:
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  // Current base point
  PointType mBasePoint;  // vector of size d = dimensionality
  // Data accumulated to represent moments
  Scalar    mN;          // 0th moment (the # of points included so far)
  PointType mDxCum;      // 1st moment (size d), the cumulative displacement
  Scalar    mDx2Cum;     // 2nd moment, the summed-square displacement
  // Storage for the inverse covariance
  Scalar    mICov;
  int       mDof;
  bool      mIsValid;
  // Temporaries (so we only allocate once)
  PointType mDx;

  /** Prepare the inverse covariance matrix
   * Taking the inverse can be subtle, particularly with regards to degenerate covariance matrices or cases where lambdaRatioThresh is non-zero.  This function implements the needed preparations.  It is called, if necessary, by any operations that need the inverse (e.g., dotProduct or statistics), so you should not generally need to call this explicitly.
   */
  void prepareInverse() {
    // Underflow issues mean that its better to compute mICov and test
    // for finiteness than it is to test whether it is nonzero and
    // then invert
    // (doing it this way fixed an extremely rare bug)
    mDof = mBasePoint.size();
    mICov = Scalar(mDof)/(mDx2Cum - mDxCum.squaredNorm()/mN);
    if (std::isfinite(mICov) && mICov > 0)
      mDof = mDxCum.size();
    else {
      mDof = 0;
      mICov = 0;
    }
    mIsValid = true;
  }
};





// Make sure Doxygen skips all the template specializations
#ifndef DOXYGEN_SHOULD_SKIP_THIS

//
// Template specialization for Diagonal covariance
//
template <typename Derived>
class AccumulatorMoments<Derived,Moments::Diagonal>
{
public:
// Typedefs and enums
  // Representation of a (vector) data point
  typedef Derived PointType;

  // The scalar, i.e., the data type for individual coordinates
  typedef typename Derived::Scalar Scalar;

  // A dynamically-sized vector, even if the class is defined using fixed size types
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> DynamicPointType;

  // A column vector of integers, used to specify groups
  typedef Eigen::Matrix<int,Derived::RowsAtCompileTime,1> ColumnOfInts;

  // The return type of the expression for the mean
  typedef typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<Scalar>, const PointType, const Eigen::CwiseUnaryOp<Eigen::internal::scalar_quotient1_op<Scalar>, const PointType> > MeanXpr;


// LIFECYCLE
  AccumulatorMoments(Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) : mThresh(lambdaRatioThresh), mIsGrouping(false) { clear(); }
  AccumulatorMoments(int d,Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) : mThresh(lambdaRatioThresh), mIsGrouping(false), mBasePoint(d), mDxCum(d), mDx2Cum(d), mICov(d), mDx(d), mTmp1(d), mTmp2(d) { clear(); }  // no grouping
  AccumulatorMoments(const ColumnOfInts& group,Scalar lambdaRatioThresh = 10*std::numeric_limits<Scalar>::epsilon()) : mThresh(lambdaRatioThresh), mIsGrouping(true), mGroup(group), mMultiplicity(group.maxCoeff()+1), mBasePoint(group.size()), mDxCum(group.size()), mDx2Cum(group.size()), mICov(group.size()), mDx(group.size()), mTmp1(group.size()), mTmp2(group.size()), mTmpG(group.maxCoeff()+1) { calculateMultiplicity(); clear(); }
   //AccumulatorMoments(const AccumulatorMoments& from);
   //~AccumulatorMoments(void);

// OPERATORS
   //AccumulatorMoments&                     operator=(const AccumulatorMoments& from);

// OPERATIONS                       
  void clear() {
    mN = 0;
    mDxCum = PointType::Zero(mDxCum.size());
    mDx2Cum = PointType::Zero(mDx2Cum.size());
    mInvState = 0;
  }

  template <typename DerivedOther>
  void restart(const Eigen::MatrixBase<DerivedOther>& coords) {
    clear();
    mBasePoint = coords;
  }

  template <typename DerivedOther>
  void addPoint(const Eigen::MatrixBase<DerivedOther>& p) {
    mN++;
    mDx = p-mBasePoint;
    mDxCum += mDx;
    mDx2Cum.array() += mDx.array().square();
    mInvState = 0;
  }

  template <typename DerivedOther>
  void addPoint(const Eigen::MatrixBase<DerivedOther>& p,Scalar w) {
    mN+=w;
    mDx = p-mBasePoint;
    mDxCum += w*mDx;
    mDx2Cum.array() += w*mDx.array().square();
    mInvState = 0;
  }

  template <typename DerivedOther>
  void removePoint(const Eigen::MatrixBase<DerivedOther>& p) {
    mN--;
    mDx = p-mBasePoint;
    mDxCum -= mDx;
    mDx2Cum.array() -= mDx.array().square();
    mInvState = 0;
  }

  template <typename DerivedOther>
  void removePoint(const Eigen::MatrixBase<DerivedOther>& p,Scalar w) {
    mN -= w;
    mDx = p-mBasePoint;
    mDxCum -= w*mDx;
    mDx2Cum.array() -= w*mDx.array().square();
    mInvState = 0;
  }

  template <typename Tl,typename Tr>
  Scalar dotProduct(const Eigen::MatrixBase<Tl>& l,const Eigen::MatrixBase<Tr>& r,int inverseType) {
    if (mInvState != inverseType)
      prepareInverse(inverseType);
    return (l.array() * r.array() * mICov.array()).sum();
  }

  template <typename T>
  void solveInPlace(Eigen::MatrixBase<T>& v,int inverseType) {
    if (mInvState != inverseType)
      prepareInverse(inverseType);
    v.array() *= mICov.array();
  }

  template <typename T>
  Scalar dotProduct(const Eigen::MatrixBase<T>& v,int inverseType) {
    if (mInvState != inverseType)
      prepareInverse(inverseType);
    return (v.array().square() * mICov.array()).sum();
  }
      
// ACCESS
  Scalar n() const { return mN; }

  const PointType& basePoint() const { return mBasePoint; }

  int dof() { if (mInvState != 1) prepareInverse(1); return mDof; }

  const PointType& inverse_diagonal() const { return mICov; }

// INQUIRY
  void statisticsNeighborhood(Scalar& T2,Scalar &dof) {
    T2 = dotProduct(mDxCum,1);
    dof = mDof;
  }

  template <typename T>
  void statisticsPoint(Scalar& T2,Scalar &dof,const Eigen::MatrixBase<T>& p) {
    if (mN == 0) {
      T2 = 0;
      dof = 0;
      return;
    }
#if PTWISESTAT_MEANCENTERED
    mDx = p - mean();
    T2 = mN * dotProduct(mDx,1);
#else
    mDx = p - mBasePoint;
    T2 = mN * dotProduct(mDx,2);
#endif
    dof = mDof;
  }
  
  MeanXpr mean() const {
    return mBasePoint+(mDxCum/mN);
  }

  Scalar meanSquareDisplacement() const { return mDx2Cum.sum()/mN; }

protected:
private:
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  // Static information
  Scalar           mThresh;
  bool             mIsGrouping;
  ColumnOfInts     mGroup;
  DynamicPointType mMultiplicity;
  // Current base point
  PointType        mBasePoint;    // size d
  // Data accumulated to represent moments
  Scalar           mN;
  PointType        mDxCum;        // size d
  PointType        mDx2Cum;       // size d (sum of squares, by coordinate)
  // Storage for the inverse covariance
  PointType        mICov;         // size d
  int              mDof;          // # of degrees of freedom
  bool             mInvState;
  // Temporaries
  PointType        mDx;           // size d
  PointType        mTmp1,mTmp2;   // size d
  DynamicPointType mTmpG;         // size nGroups

  // Function to calculate sums of grouped coordinates of v
  void accumulateGroup(DynamicPointType &vg,const PointType &v) {
    vg = DynamicPointType::Zero(vg.size());
    for (int i = 0; i < mGroup.size(); i++)
      vg[mGroup[i]] += v[i];
  }

  // Function to calculate the number of coordinates in each group
  void calculateMultiplicity() {
    mMultiplicity = DynamicPointType::Zero(mMultiplicity.size());
    for (int i = 0; i < mGroup.size(); i++)
      mMultiplicity[mGroup[i]]++;
  }

  // Calculate the inverse matrix. inverseType = 1 corresponds to
  // using the sample covariance (mean-centered), inverseType = 2
  // corresponds to the basepoint-centered covariance.
  void prepareInverse(int inverseType) {
    if (inverseType == 1) {
      mTmp1.array() = mDx2Cum.array() - mDxCum.array().square()/mN;
      mTmp1.array() = (mTmp1.array() < 0).select(0,mTmp1.array());
      prepareInverseCalc(mTmp1);
    } else
      prepareInverseCalc(mDx2Cum);
    mInvState = inverseType;
  }

  // Function to calculate the inverse matrix. We pass in a vector
  // that is either mDx2Cum (for basepoint-centered covariance) or
  // mDx2Cum - mDxCum.^2/N (for mean-centered)
  void prepareInverseCalc(const PointType& diag) {
    Scalar thresh;
    int i,j;
    mDof = 0;
    if (mIsGrouping) {
      // Sum the coordinate groups
      accumulateGroup(mTmpG,diag);
      mTmpG.array() = mTmpG.array()/mMultiplicity.array();
      // Keep only those groups that exceed threshold
      thresh = mTmpG.maxCoeff() * mThresh;
      // Calculate inverse and split back out into individual coordinates
      for (i = 0; i < mICov.size(); i++) {
	j = mGroup[i];
	if (mTmpG[j] > thresh) {
	  mICov[i] = 1.0/mTmpG[j];
	  mDof++;  // keep track of dof
	}
	else
	  mICov[i] = 0;
      }
    } else {
      // Keep only coordinates that exceed threshold
      thresh = diag.maxCoeff() * mThresh;
      for (i = 0; i < mICov.size(); i++) {
	if (diag[i] > thresh) {
	  mICov[i] = 1.0/diag[i];
	  mDof++;
	}
	else
	  mICov[i] = 0;
      }
    }
  }

};







//
// Template specialization for Full covariance
//
template <typename Derived>
class AccumulatorMoments<Derived,Moments::Full>
{
public:
// Typedefs and enums
  // Representation of a (vector) data point
  typedef Derived PointType;

  // The scalar, i.e., the data type for individual coordinates
  typedef typename Derived::Scalar Scalar;

  // The type of the covariance matrix
  typedef typename Eigen::Matrix<Scalar,PointType::RowsAtCompileTime,PointType::RowsAtCompileTime> MatrixType;

  // The type of the covariance matrix decomposition
  typedef typename Eigen::LDLT<MatrixType> CovarType;

  // The return type of the expression for the mean
  typedef typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<Scalar>, const PointType, const Eigen::CwiseUnaryOp<Eigen::internal::scalar_quotient1_op<Scalar>, const PointType> > MeanXpr;


// LIFECYCLE
  AccumulatorMoments(Scalar lambdaRatioThresh = 1000*std::numeric_limits<Scalar>::epsilon()) : mThresh(lambdaRatioThresh) { clear(); }
  AccumulatorMoments(int d,Scalar lambdaRatioThresh = 1000*std::numeric_limits<Scalar>::epsilon()) : mThresh(lambdaRatioThresh), mBasePoint(d), mDxCum(d), mDx2Cum(d), mCov(d), mDx(d), mID(d), mTmp1(d), mTmp2(d) { clear(); }
  //AccumulatorMoments(const ColumnOfInts& group,Scalar lambdaRatioThresh = 0) : mThresh(lambdaRatioThresh), mIsGrouping(true), mGroup(group), mMultiplicity(group.maxCoeff()), mBasePoint(group.size()), mDxCum(group.size()), mDx2Cum(group.size()), mDx(group.size()), mTmpG1(group.maxCoeff()), mTmpG2(group.maxCoeff()) { calculateMultiplicity(); clear(); }
   //AccumulatorMoments(const AccumulatorMoments& from);
   //~AccumulatorMoments(void);

// OPERATORS
   //AccumulatorMoments&                     operator=(const AccumulatorMoments& from);  

// OPERATIONS                       
  void clear() {
    mN = 0;
    mDxCum = PointType::Zero(mDxCum.size());
    mDx2Cum.clear();
    mIsValid = false;
  }

  template <typename DerivedOther>
  void restart(const Eigen::MatrixBase<DerivedOther>& coords) {
    clear();
    mBasePoint = coords;
  }

  template <typename DerivedOther>
  void addPoint(const Eigen::MatrixBase<DerivedOther>& p) {
    mN++;
    mDx = p-mBasePoint;
    mDxCum += mDx;
    mDx2Cum.rankUpdate(mDx);
    mIsValid = false;
  }

  template <typename DerivedOther>
  void addPoint(const Eigen::MatrixBase<DerivedOther>& p,Scalar w) {
    mN += w;
    mDx = p-mBasePoint;
    mDxCum += w*mDx;
    mDx2Cum.rankUpdate(mDx,w);
    mIsValid = false;
  }

  template <typename DerivedOther>
  void removePoint(const Eigen::MatrixBase<DerivedOther>& p) {
    mN--;
    mDx = p-mBasePoint;
    mDxCum -= mDx;
    mDx2Cum.rankUpdate(mDx,-1);
    mIsValid = false;
  }

  template <typename DerivedOther>
  void removePoint(const Eigen::MatrixBase<DerivedOther>& p,Scalar w) {
    mN -= w;
    mDx = p-mBasePoint;
    mDxCum -= w*mDx;
    mDx2Cum.rankUpdate(mDx,-w);
    mIsValid = false;
  }

  template <typename Tl,typename Tr>
  Scalar dotProduct(const Eigen::MatrixBase<Tl>& l,const Eigen::MatrixBase<Tr>& r) {
    if (!mIsValid)
      prepareInverse();
    mTmp1 = mCov.transpositionsP() * l;
    mCov.matrixL().solveInPlace(mTmp1);
    mTmp2 = mCov.transpositionsP() * r;
    mCov.matrixL().solveInPlace(mTmp2);
    return (mTmp1.array() * mTmp2.array() * mID.array()).sum();
  }

  template <typename T>
  void solveInPlace(Eigen::MatrixBase<T>& v) {
    if (!mIsValid)
      prepareInverse();
    mTmp1 = mCov.transpositionsP()*v;
    mCov.matrixL().solveInPlace(mTmp1);
    v.array() = mTmp1.array() * mID.array();
    mCov.matrixU().solveInPlace(v);
  }

  template <typename T>
  Scalar dotProduct(const Eigen::MatrixBase<T>& v) {
    if (!mIsValid)
      prepareInverse();
    mTmp1 = mCov.transpositionsP() * v;
    mCov.matrixL().solveInPlace(mTmp1);
    return (mTmp1.array().square() * mID.array()).sum();
  }
  
// ACCESS
  Scalar n() const { return mN; }

  const PointType& basePoint() const { return mBasePoint; }

  int dof() { if (!mIsValid) prepareInverse(); return mDof; }

  const CovarType& covariance() const { return mCov; }

  const PointType& inverse_diagonal() const { return mID; }

// INQUIRY
  void statisticsNeighborhood(Scalar& T2,Scalar &dof) {
    T2 = dotProduct(mDxCum);
    dof = mDof;
  }

  template <typename T>
  void statisticsPoint(Scalar& T2,Scalar &dof,const Eigen::MatrixBase<T>& p) {
    if (mN == 0) {
      T2 = 0;
      dof = 0;
      return;
    }
#if PTWISESTAT_MEANCENTERED
    mDx = p - mean();
    T2 = mN * dotProduct(mDx);
#else
    // Calculate (p-basepoint)^T (C + X X^T/N)^{-1} (p-basepoint) by Sherman-Morrison
    mDx = p - mBasePoint;
    T2 = mN * dotProduct(mDx);
    Scalar T2X = dotProduct(mDxCum);
    Scalar dp = dotProduct(mDx,mDxCum);
    T2 -= dp*dp/(1+T2X/mN);
#endif
    dof = mDof;
  }
  
  MeanXpr mean() const {
    return mBasePoint+(mDxCum/mN);
  }

  Scalar meanSquareDisplacement() const { return mDx2Cum.vectorD().sum()/mN; }


protected:
private:
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  // Static information
  Scalar           mThresh;
  // Current base point
  PointType        mBasePoint;    // size d
  // Data accumulated to represent moments
  Scalar           mN;
  PointType        mDxCum;        // size d
  CovarType        mDx2Cum;       // size dxd (sometimes called the SSP, sum of squares and products, matrix)
  CovarType        mCov;
  // Storage for the inverse covariance
  PointType        mID;           // size d (inverse of diagonals)
  int              mDof;          // # of degrees of freedom
  bool             mIsValid;
  // Temporaries
  PointType        mDx;           // size d
  PointType        mTmp1,mTmp2;   // size d

  // Note: we always invert the sample covariance, even when
  // PTWISESTAT_MEANCENTERED is true. We always want the sample
  // covariance for the neighborhood statistics, and it is not
  // numerically stable to compute the neighborhood statistics using
  // Sherman-Morrison. So, when PTWISESTAT_MEANCENTERED is true, we
  // compute the basepoint-centered dot product using traditional
  // Sherman-Morrison. When the sample covariance is rank-deficient,
  // this is not accurate in a linear-algebra sense, but we could make
  // it so using a modified formula (see Riedel 1991, A Sherman
  // Morrison Woodbury identity for rank augmenting matrices with
  // application to centering).  However, we don't use this modified
  // formula, because in this case I think the traditional
  // Sherman-Morrison better corresponds to what we want to do: to the
  // extent that the difference between the basepoint and the mean has
  // not been "fleshed out" in the sample covariance, just ignore this
  // difference.  Since doing a pseudoinverse on the sample covariance
  // throws out any undetermined component of the difference between
  // the basepoint and the mean, this better achieves our goal.
  void prepareInverse() {
    Scalar thresh,tmp;
    int i;
    mDof = 0;
    mCov = mDx2Cum;
    mCov.rankUpdate(mDxCum,Scalar(-1)/mN);
    thresh = mCov.vectorD().maxCoeff() * mThresh;
    thresh = std::max(Scalar(0),thresh);
    for (i = 0; i < mDx.size(); i++) {
      tmp = mCov.vectorD().coeff(i);
      if (tmp > thresh) {
	mID[i] = 1.0/tmp;
	mDof++;
      }
      else
	mID[i] = 0;
    }
    mIsValid = true;
  }
};



#endif // DOXYGEN_SHOULD_SKIP_THIS  (for template specializations)


#endif  // AccumulatorMoments_h
