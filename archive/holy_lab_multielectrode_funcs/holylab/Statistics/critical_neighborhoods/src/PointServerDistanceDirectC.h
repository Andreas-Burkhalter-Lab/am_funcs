#ifndef PointServerDistanceDirectC_h
#define PointServerDistanceDirectC_h

// SYSTEM INCLUDES
//
#include <vector>
#include <algorithm>
#include "Eigen/Dense"

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "Metric.h"
#include "IndirectCompare.h"
#include "OutlierStatistics.h"

// FORWARD REFERENCES
//

/**  \brief A simple PointServer ordering points by distance
 *
 * \details #include "PointServerDistanceDirectC.h" <BR>
 *
 * This serves points in order of increasing distance from a base point.  Every point is checked, so it is an \f$O(N)\f$ algorithm.
 * \tparam Derived The matrix (or array) used to represent the collection of data points (on the columns of the matrix)
 * \tparam MetricType A class supporting a general interface for measuring distances between points.  See Metric.  The default metric is Euclidean.
 *  
 */

template <typename Derived,typename MetricType>
class PointServerDistanceDirectC
{
public:
  /** The scalar data type
   */
  typedef typename Derived::Scalar Scalar;
  /** The type of a single data point (a column), with allocated storage
   */
  typedef typename Eigen::Matrix<Scalar, Derived::RowsAtCompileTime, 1> PointType;
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
#endif // DOXYGEN_SHOULD_SKIP_THIS

// LIFECYCLE
  /** The constructor.
   * \note After construction, you should use basePoint(bp) to set or change the base point; the PointServer is not ready to use until you set the base point (construction does not compute distances to all points, but setting the base point does).
   * \param X A matrix with each data point corresponding to a column of X
   * \param metric An instantiation of the metric, e.g., Metric::Euclidean<PointType> metric. Be aware that any parameters of the metric may change as a consequence to the call to update.
   */
  PointServerDistanceDirectC(const Derived &X, MetricType &metric, OutlierStatistics<Scalar>& os) : mrX(X), mSortOrder(X.cols()), mDistance(X.cols()), mBasePoint(X.col(0)), mIndex(0), mrMetric(metric), mCheckOutlier(true), mpOutlier(&os) { ; }
  PointServerDistanceDirectC(const Derived &X, MetricType &metric) : mrX(X), mSortOrder(X.cols()), mDistance(X.cols()), mBasePoint(X.col(0)), mIndex(0), mrMetric(metric), mCheckOutlier(false), mpOutlier(NULL) { ; }
  // Prevent copying: private copy constructor & assignment operator
  // Default destructor
// OPERATORS
  /** Prefix increment operator
   * \return Nothing
   */
  void operator++() { mIndex++; }

  /** Postfix increment operator
   * \return Nothing
   */
  void operator++(int unused) { mIndex++; }
  

// OPERATIONS
  /** Restart at the beginning, keeping the current base point.  Consequently, this has no computational overhead.
   */
  void restart() { mIndex = 0; }

  /** Set the base point.
   * This triggers computation of the distance between the new base point and all points in the data set, and then ordering by distance.  Consequently, it can be computationally intensive.
   * \param bp A reference to the coordinates of the new base point
   */
  template <typename T> void basePoint(const Eigen::MatrixBase<T>& bp) { mBasePoint = bp; orderPoints();}

  /** Set the base point to one of the data points.
   * This triggers computation of the distance between the new base point and all points in the data set, and then ordering by distance.  Consequently, it can be computationally intensive.
   * \param bp A reference to the coordinates of the new base point
   */
  void basePointIndex(int index) { mBasePoint = point(index); orderPoints();}

  template <typename Accumulator>
  void update(const Accumulator& acc) {
    basePoint(acc.mean());
    mrMetric.update(acc);
  }
    
  
// ACCESS
  /** Obtain the number of points in the sequence 
   * \return The number of points in the data set
   */
  int N() const { return mSortOrder.size(); }

  /** Obtain the number of points.
   * \return The number of points in the data set
   */
  int Ntot() const { return mrX.cols(); }

  /** Obtain the current point's index.
   * \return The integer index of the current point (zero-offset)
   */
  int currentIndex() const { return mSortOrder[mIndex]; }

  /** Retrieve the current data point.
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr currentPoint() const { return mrX.col(currentIndex()); }

  /** Retrieve a particular data point.
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr point(int index) const { return mrX.col(index); }

  /** Retrieve the base point.
   * \return A constant reference to the coordinates of the base point
   */
  const PointType& basePoint() const { return mBasePoint; }


// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence.
   */
  bool isAtEnd() const { return (mIndex >= mrX.cols()); }

  /** Test whether the current point "ties" the previous one, meaning it has the same distance.
   */
  bool isTiesPrevious() const {
    if (mIndex > 0)
      return (mDistance[mSortOrder[mIndex]] == mDistance[mSortOrder[mIndex-1]]);
    else
      return false;
  }

protected:
private:
  const Derived       &mrX;      // a reference to the data points
  std::vector<int>    mSortOrder;// the sequence vector
  std::vector<Scalar> mDistance; // measure of distance
  PointType           mBasePoint;// the base point from which displacements are measured
  int                 mIndex;    // the current position in the sequence vector
  MetricType          &mrMetric; // the metric to measure distance
  bool                mCheckOutlier; // true if we should check for outliers
  OutlierStatistics<Scalar>* mpOutlier;// reference to class for outlier detection

  void orderPoints() {
    // Calculate the distance from the base point to all data points
    //Eigen::Map<const Derived> mX(pX,d,N);
    int i;
    if (mCheckOutlier) {
      // If we are doing outlier detection, we have to use the
      // squared-distance, otherwise the statistics will be wrong
      for (i = 0; i < mrX.cols(); i++) {
	mDistance[i] = mrMetric.squaredDistance(mBasePoint,mrX.col(i));
	mSortOrder[i] = i;
      }
    }
    else {
      // If we're not doing outlier detection, just use the fastest
      // metric available
      for (i = 0; i < mrX.cols(); i++) {
	mDistance[i] = mrMetric.fastDistance(mBasePoint,mrX.col(i));
	mSortOrder[i] = i;
      }
    }
    // Sort in increasing order of distance
    IndirectLess<std::vector<Scalar> > comp(mDistance);
    std::sort(mSortOrder.begin(),mSortOrder.end(),comp);
    // Look for outliers
    if (mCheckOutlier) {
      int n = mpOutlier->checkSequence(mSortOrder,mDistance,mrX.rows());
      mSortOrder.resize(n);
    }
    restart();
  }

  // Block copying
  PointServerDistanceDirectC(const PointServerDistanceDirectC& from);
  PointServerDistanceDirectC&                     operator=(const PointServerDistanceDirectC& from);
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerDistanceDirectC_h

