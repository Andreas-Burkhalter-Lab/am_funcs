#ifndef PointServerFilter1d_h
#define PointServerFilter1d_h

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
#include "OutlierStatistics.h"

// FORWARD REFERENCES
//

/**  \brief A simple PointServer ordering points by distance
 *
 * \details #include "PointServerFilter1d.h" <BR>
 *
 * This serves points in order of increasing time separation from a particular moment in time.  It can be "causal" (the default) or "acausal."  With "causal," only points from an earlier time are served.  With "acausal," points are served symmetrically around the current time; first the left point will be served, then the right point will be served (increment and check the isTiesPrevious() method to determine whether one should consider the points on equal footing).
 * \tparam Derived The vector (or matrix) used to represent the value(s) measured at all time points (if multiply-valued, each time point is a column of this matrix).
 * \tparam acausal (default false) A boolean; if true, the points are served symmetrically around the "base time."
 * \tparam MetricType A class supporting a general interface for measuring distances between points (see Metric).  This is used for outlier detection (if desired), not for sequencing the points.  Default is Euclidean.
 */

template <typename Derived,bool acausal = false,typename MetricType = Metric::Euclidean>
class PointServerFilter1d
{
public:
  /** The scalar data type  */
  typedef typename Derived::Scalar Scalar;
  /** The type of a single data point (a column), with allocated storage
   */
  typedef typename Eigen::Matrix<Scalar, Derived::RowsAtCompileTime, 1> PointType;
  /** The type of the input data */
  typedef Derived MatrixType;
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
#endif // DOXYGEN_SHOULD_SKIP_THIS

// LIFECYCLE
  /** The constructor.
   * \note After construction, you should use basePointIndex(index) to set or change the base point; the PointServer is not ready to use until you set the base point.
   * \param X A vector or matrix, with each "time point" corresponding to a column of X
   * \param os An instantiation of class OutlierStatistics, supply this only if you want to do outlier detection.
   * \param metric An instantiation of the metric, e.g., Metric::Euclidean<PointType> metric. Use only for outlier detection.  Be aware that any parameters of the metric may change as a consequence to the call to update.
   */
  PointServerFilter1d(const Derived &X, OutlierStatistics<Scalar>& os, MetricType &metric) : mrX(X), mBasePoint(X.col(0)), mIndex(0), mCheckOutlier(true), mpOutlier(&os), mpMetric(&metric) { ; }
  PointServerFilter1d(const Derived &X) : mrX(X), mBasePoint(X.col(0)), mIndex(0), mCheckOutlier(false), mpOutlier(NULL), mpMetric(NULL) { ; }
  // Prevent copying: private copy constructor & assignment operator
  // Default destructor
// OPERATORS
  /** Prefix increment operator
   * \return Nothing
   */
  void operator++() { inc(); }

  /** Postfix increment operator
   * \return Nothing
   */
  void operator++(int unused) { inc(); }
  

// OPERATIONS
  /** Restart at the beginning, keeping the current base point.
   */
  void restart() {
    mIndex = mBasePointIndex;
    mN = 0;
    // Tasks involved in outlier detection
    mOutlierTriggered = false;
    if (mCheckOutlier) {
      mR2cum = mpMetric->squaredDistance(mBasePoint,mrX.col(mBasePointIndex));
    }
  }

  /** Set the base point value.
   * This sets the value associated with the current base point; it does not change the "time" at which the base point is defined, and hence does not change the ordering of points.
   * \param bp A reference to the coordinates of the new base point
   */
  template <typename T> void basePoint(const Eigen::MatrixBase<T>& bp) { mBasePoint = bp; restart(); }

  /** Set the "time point" of the base point to one of the data points.
   * \param bp A reference to the coordinates of the new base point
   * \note The "value" of the base point is set equal to the value of the indexed point.
   */
  void basePointIndex(int index) { mBasePointIndex = index; mBasePoint = point(index); restart(); }

  /** Reset the accumulator
   * \param acc The accumulator
   */
  template <typename AccumulatorType>
  void update(const AccumulatorType& acc) {
    //std::cout << acc.n() << " points, mean " << acc.mean() << std::endl;
    basePoint(acc.mean());
    mpMetric->update(acc);
  }
    
  /** Add the current point to the accumulator
   * \param acc the accumulator
   */
  template <typename AccumulatorType>
  void addCurrentPoint(AccumulatorType& acc) {
    acc.addPoint(currentPoint());
  }

  
// ACCESS
  /** Obtain the number of points in the sequence.
   * \return The number of points in the sequence
   */
  int N() const {
    if (acausal)
      return mrX.cols();
    else
      return mBasePointIndex+1;
  }

  /** Obtain the number of points in the whole data set.
   * \return The number of points in the data set
   */
  int Ntot() const { return mrX.cols(); }

  /** Obtain the current point's index.
   * \return The integer index of the current point (zero-offset)
   */
  int currentIndex() const { return mIndex; }

  /** Retrieve the current value(s).
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr currentPoint() const { return mrX.col(currentIndex()); }

  /** Retrieve the value(s) at a particular time.
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr point(int index) const { return mrX.col(index); }

  /** Retrieve the base point.
   * \return A constant reference to the coordinates of the base point
   */
  const PointType& basePoint() const { return mBasePoint; }

  /** Retrieve the base point's "time".
   * \return The integer index of the "time" for the current base point
   */
  int basePointIndex() const { return mBasePointIndex; }

// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence.
   */
  bool isAtEnd() const { return (mIndex < 0 || mIndex >= mrX.cols() || mOutlierTriggered); }

  /** Test whether the current point "ties" the previous one, meaning it has the same distance.
   */
  bool isTiesPrevious() const {
    return (acausal && mIndex > mBasePointIndex && mIndex < 2*mBasePointIndex);
  }

  /** Filtering point servers do not accept weights */
  bool hasWeight() const { return false; }
  int weight(int index) const { return 1; }

protected:
private:
  const Derived &mrX;      // a reference to the data points
  int           mBasePointIndex; // base point "time" (determines ordering)
  PointType     mBasePoint;// base point "value"
  int           mIndex;    // the current position in the sequence vector
  int           mN;         // the number of points included so far
  bool          mCheckOutlier; // true if we should check for outliers
  OutlierStatistics<Scalar> *mpOutlier; // for outlier detection
  bool          mOutlierTriggered;
  Scalar        mR2cum;
  MetricType    *mpMetric; // the metric to measure distance
  
  void inc() {
    int i;
    if (acausal) {
      // We're alternating between earlier & later neighbors
      if (mIndex < mBasePointIndex) {
	// We just did the earlier one, now check the later one
	i = 2*mBasePointIndex-mIndex;
	if (i < mrX.cols())
	  mIndex = i;  // it was within bounds, keep it
	else
	  mIndex--;  // later was out of bounds, use a more distant earlier nbr
      } else {
	// We just did the later one, switch to a more-distant earlier one
	i = 2*mBasePointIndex-mIndex-1;
	if (i >= 0)
	  mIndex = i;  // it was within bounds, keep it
	else
	  mIndex++; // earlier was out of bounds, use a more distant later nbr
      }
    } else
      mIndex--;  // causal case (go to an earlier time)
    if (mCheckOutlier) {
      Scalar r2 = mpMetric->squaredDistance(mBasePoint,mrX.col(mIndex));
      mOutlierTriggered = mOutlierTriggered | mpOutlier->isSignificant(r2/(mR2cum/mN),1,mN);
      mR2cum += r2;
    }
    mN++;
  }
  

  // Prevent copying
  PointServerFilter1d(const PointServerFilter1d& from);
  PointServerFilter1d& operator=(const PointServerFilter1d& from);
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerFilter1d_h

