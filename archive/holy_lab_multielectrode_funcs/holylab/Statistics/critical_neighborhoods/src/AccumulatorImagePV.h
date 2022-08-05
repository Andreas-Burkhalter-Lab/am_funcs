/** \brief  Accumulate statistics of pixel spatial and value displacements
 *
 * #include "AccumulatorImagePV.h" <BR>
 *
 * \tparam CoordP The vector type used to represent the spatial position.  This must have its own allocated storage (e.g, Eigen::Matrix). This should be a floating-point type (float/double).
 * \tparam CoordV The value type holding the image information at each pixel. For example, for a 16-bit grayscale image this would be Eigen::Matrix<uint16_t,1,1>, for an RGB image it might be a Vector3f, etc.  It is OK to use an integer type here; internally, integer types will be cast to the same scalar type as used for CoordP.
 * \tparam ctype The enum type used to indicate the covariance model for the value data (see namespace Moments)
 */
 
#ifndef AccumulatorImagePV_h
#define AccumulatorImagePV_h

// SYSTEM INCLUDES
//
#include "Eigen/Dense"

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "AccumulatorMoments.h"
#include "PointServerImagePV.h"

// FORWARD REFERENCES
//



template <typename CoordP,typename CoordV,Moments::CovarianceModel ctype = Moments::Isotropic>
class AccumulatorImagePV
{
public:
  /** Typedef for the scalar, i.e., the data type for individual coordinates */
  typedef typename CoordP::Scalar Scalar;

private:
  typedef AccumulatorMoments<CoordP,Moments::Isotropic>     AccPType;
  // Value coordinates, cast to type Scalar
  typedef Eigen::Matrix<Scalar,CoordV::RowsAtCompileTime,1> CoordVP;
  typedef AccumulatorMoments<CoordVP,ctype>                 AccVType;
  // Return type of mean
  typedef ImagePoint<const typename AccPType::MeanXpr,const typename AccVType::MeanXpr> MeanXpr;
public:
// LIFECYCLE
  /** Default constructor
   *  This is the "typical" constructor when the dimensionality is specified at compile time.
   *  \param lambdaRatioThresh Any eigenvalues/eigenvectors of the "value" covariance matrix satisfying \f$\lambda_i/\lambda_\mathrm{max} < \f$ lambdaRatioThresh will be excluded (default value 0).
   */
  AccumulatorImagePV(Scalar lambdaRatioThresh = 0) : mAccP(), mAccV(lambdaRatioThresh) { clear(); }
  /** Constructor specifying dimensionality
   *  This is the "typical" constructor when the dimensionality is specified at run time.
   * \param dSpatial The spatial dimensionality
   * \param dValue The number of values used to represent each pixel (e.g., 1 for grayscale and 3 for RGB)
   *  \param lambdaRatioThresh Any eigenvalues/eigenvectors of the "value" covariance matrix satisfying \f$\lambda_i/\lambda_\mathrm{max} < \f$ lambdaRatioThresh will be excluded (default value 0).
   */
  AccumulatorImagePV(int dSpatial,int dValue, Scalar lambdaRatioThresh = 0) : mAccP(dSpatial), mAccV(dValue,lambdaRatioThresh), mTmpVP(dValue) { clear(); }


// OPERATORS

// OPERATIONS                       
  /** Restart from the current basepoint
   */
  void clear() { mAccP.clear(); mAccV.clear(); }

  /** Restart with a new basepoint
   *
   * \param bp The coordinates of the new base point
   */
  template <typename PT,typename VT> void restart(const ImagePointRef<PT,VT>& bp) {
    //std::cout << "restart:\n  position " << bp.position().transpose() << "\n  value " << bp.value().transpose() << std::endl;
    mAccP.restart(bp.position());
    mAccV.restart(bp.value());
  }

  /** Add a new point to the neighborhood, and incorporate into moments
   *
   * \param p A reference to the point
   */
  template <typename PT,typename VT> void addPoint(const ImagePoint<PT,VT>& p) {
    mAccP.addPoint(p.position());
    // Must cast to floating-point type so don't overflow
    mTmpVP = p.value().template cast<Scalar>();
    mAccV.addPoint(mTmpVP);
    //std::cout << "added position " << p.position().transpose() << ", value " << mTmpVP.transpose() << std::endl;
  }

  /** Remove a point from the neighborhood, subtracting its effect on moments
   *
   * \param p A reference to the point
   */
  template <typename PT,typename VT> void removePoint(const ImagePoint<PT,VT>& p) {
    mAccP.removePoint(p.position());
    // Must cast to floating-point type so don't overflow
    mTmpVP = p.value().template cast<Scalar>();
    mAccV.removePoint(mTmpVP);
    //std::cout << "removed position " << p.position().transpose() << ", value " << mTmpVP.transpose() << std::endl;
  }

  template <typename PT,typename VT> void removePoint(const ImagePointRef<PT,VT>& p) {
    mAccP.removePoint(p.position());
    // Must cast to floating-point type so don't overflow
    mTmpVP = p.value().template cast<Scalar>();
    mAccV.removePoint(mTmpVP);
  }

  /** In preparation for statistics */
  void prepareInverse() {
    mAccP.prepareInverse();
    mAccV.prepareInverse();
  }


// ACCESS
  /** 
   *
   * \return The number of points accumulated so far
   */
  int n() const { return mAccP.n(); }

  /** Access the position accumulator
   *
   * \return A constant reference to the position Accumulator
   */
  const AccPType& position() const { return mAccP; }

  /** Access the value accumulator
   *
   * \return A constant reference to the value Accumulator
   */
  const AccVType& value() const { return mAccV; }

// INQUIRY
  /** Calculate the terms needed to determine statistical signficance
   *
   * \param chisq A reference for the storage of \f$\chi^2\f$
   * \param dof A reference for the storage of the number of degrees of freedom
   */
  void statistics(Scalar& chisqNeighborhood,Scalar& chisqPoint,Scalar &dof) {
    Scalar chisqN2,chisqP2,dof2;
    mAccP.statistics(chisqNeighborhood,chisqPoint,dof);
    mAccV.statistics(chisqN2,chisqP2,dof2);
    chisqNeighborhood += chisqN2;
    chisqPoint += chisqP2;
    dof += dof2;
  }

  void statistics(Scalar& chisqNeighborhood,Scalar &dof) {
    Scalar chisqN2,dof2;
    mAccP.statistics(chisqNeighborhood,dof);
    mAccV.statistics(chisqN2,dof2);
    chisqNeighborhood += chisqN2;
    dof += dof2;
  }

  template <typename PT,typename VT>
  Scalar statistics(const ImagePoint<PT,VT>& ip) {
    return mAccP.statistics(ip.position()) + mAccV.statistics(ip.value().template cast<Scalar>());
  }

  template <typename PT,typename VT>
  Scalar statistics(const ImagePointRef<PT,VT>& ip) {
    return mAccP.statistics(ip.position()) + mAccV.statistics(ip.value().template cast<Scalar>());
  }

  /** Calculate the center of mass of all points accumulated thus far
   *
   * \return An ImagePoint expressing the position and value
   */
  MeanXpr mean() const { return MeanXpr(mAccP.mean(),mAccV.mean()); }

  /** Calculate the ratio of variances
   *
   * \return A scalar equal to PositionMSD/ValueMSD (MSD = mean square displacement)
   */
  Scalar varianceRatio() const { return mAccP.meanSquareDisplacement()/mAccV.meanSquareDisplacement(); }

  void print() const {
    std::cout << "position: ";
    mAccP.print();
    std::cout << "value: ";
    mAccV.print();
  }

protected:
private:
  AccPType mAccP;
  AccVType mAccV;
  enum { NeedsToAlign = (sizeof(CoordVP)%16)==0};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
  CoordVP  mTmpVP;
};


#endif  // AccumulatorImagePV_h
