/**  \brief A PointServer in which order is specified "by hand".
 *
 * \details #include "PointServerPreorderedC.h" <BR>
 *
 * This is perhaps the simplest PointServer, for which you to specify the explicit order in which points are served.
 * \tparam Derived The matrix (or array) used to represent the collection of data points (on the columns of the matrix)
 * \tparam OrderType A vector type which specifies the sequence as the column index of the data matrix. The only operations that the vector sortOrder needs to support are operator[] (as in sortOrder[i]) and size().
 * \tparam offset An optional integer, use 1 if the values in sortOrder are unit-offset (default 0 = zero-offset). Having this be a template parameter should mean that no overhead is incurred when this is set to 0.
 *  
 */

#ifndef PointServerPreorderedC_h
#define PointServerPreorderedC_h

// SYSTEM INCLUDES
//
#include "Eigen/Dense"

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//


template <typename Derived,typename OrderType,int offset = 0>
class PointServerPreorderedC
{
public:
// Typedefs and enums
  /** The scalar data type */
  typedef typename Derived::Scalar Scalar;
      
  /** The basic type used to represent a single point (a column-vector) */
  typedef typename Eigen::Matrix<Scalar,Derived::RowsAtCompileTime,1> PointType;

  /** Flag to indicate whether weights are supplied (available at compile time) */
  enum { isWeighted = 0 };

// LIFECYCLE
  /** The constructor.
   * By default the base point is set to the first data point, as indicated by sortOrder[0].  You can use basePoint(bp) to change the base point.
   * \param X A matrix with each data point corresponding to a column of X
   * \param sortOrder A vector listing the sequence of points, expressed as column numbers of X. If desired, this vector can be shorter than the number of data points.
   */
  PointServerPreorderedC(const Derived &X,const OrderType &sortOrder) : mrX(X), mrSortOrder(sortOrder), mBasePoint(X.col(sortOrder[0]-offset)), mIndex(0) {;}

  // Default destructor
  //~PointServerPreorderedC;

  // Prevent copying: copy constructor and assignment operator are private

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
  /** Restart at the beginning
   */
  void restart() { mIndex = 0; }

  /** Add point to accumulator
   * \param index the index of the point to add
   * \param acc the accumulator
   */
  template <typename AccumulatorType>
  void addPoint(int index,AccumulatorType& acc) const { acc.addPoint(point(index)); }
  
  /** Add the current point to accumulator
   * \param acc the accumulator
   */
  template <typename AccumulatorType>
  void addCurrentPoint(AccumulatorType& acc) const { acc.addPoint(currentPoint()); }

// ACCESS
  /** Obtain the total number of points in the sequence
   * \return An integer with the total number of points to be considered
   */
  int N() const { return mrSortOrder.size(); }

  /** Obtain the total number of data points
   * \return An integer with the total number of data points
   */
  int Ntot() const { return mrX.cols(); }

  /** Obtain the current point's index
   * \return The integer index of the current point (zero-offset)
   */
  int currentIndex() const { return mrSortOrder[mIndex] - offset; }

  /** Retrieve the current point
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr currentPoint() const { return mrX.col(currentIndex()); }

  /** Retrieve an arbitrary point
   * \return A constant reference to the coordinates of the current point
   */
  typename Derived::ConstColXpr point(int index) const { 
    //mexPrintf("Fetching column %d\n",index);
    return mrX.col(index);
  }

  /** "Retrieve" the weight (always 1 for the unweighted class)
   * \return The weight
   */
  Scalar weight(int index) const { return 1; }

  /** Retrieve the base point.
   * Note that while the base point does not affect the order of points, it is needed, e.g., to set up the Accumulator in CriticalNeighborhood::triggerCriterion.
   * \return A constant reference to the coordinates of the base point
   */
  const PointType& basePoint() const { return mBasePoint; }

  /** Set the base point
   * \param bp A reference to the coordinates of the new base point
   */
  template <typename T> void basePoint(const Eigen::DenseBase<T>& bp) { mBasePoint = bp; }


// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence
   */
  bool isAtEnd() const { return (mIndex >= mrSortOrder.size()); }

  /** Test whether the current point "ties" the previous one. Always false.
   */
  bool isTiesPrevious() const { return false; }


protected:
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  enum { NeedsToAlign = (sizeof(PointType)%16)==0 };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
#endif // DOXYGEN_SHOULD_SKIP_THIS

  const Derived&   mrX;          // a reference to the data matrix
  const OrderType& mrSortOrder;  // a reference to the sequence vector
  PointType        mBasePoint;   // from which displacements are measured
  int              mIndex;       // the current position in the sequence vector

private:
  // Make these private until we know we need them
  PointServerPreorderedC(const PointServerPreorderedC& from);
  PointServerPreorderedC& operator=(const PointServerPreorderedC& from);
};







/** A variant in which each data point is given a scalar weight */
template <typename Derived,typename DerivedWeight,typename OrderType,int offset = 0>
class PointServerPreorderedCWeighted : public PointServerPreorderedC<Derived,OrderType,offset> {
public:
// TYPEDEFS AND ENUMS
  typedef PointServerPreorderedC<Derived,OrderType,offset> BaseType;
  /** The scalar data type */
  typedef typename Derived::Scalar Scalar;
  /** Flag to indicate whether weights are supplied (available at compile time) */
  enum { isWeighted = 1 };

// LIFECYCLE
  /** The constructor.
   * By default the base point is set to the first data point, as indicated by sortOrder[0].  You can use basePoint(bp) to change the base point.
   * \param X A matrix with each data point corresponding to a column of X
   * \param w A vector of scalar weights, one per data point
   * \param sortOrder A vector listing the sequence of points, expressed as column numbers of X. If desired, this vector can be shorter than the number of data points.
   */
  PointServerPreorderedCWeighted(const Derived &X,const DerivedWeight &w,const OrderType &sortOrder) : BaseType(X,sortOrder), mrWeight(w) {;}

// OPERATIONS
  /** Add weighted point to accumulator
   * \param index the index of the point to add
   * \param acc the accumulator
   */
  template <typename AccumulatorType>
  void addPoint(int index,AccumulatorType& acc) const { acc.addPoint(this->point(index),weight(index)); }
  
  /** Add the current point to accumulator
   * \param acc the accumulator
   */
  template <typename AccumulatorType>
  void addCurrentPoint(AccumulatorType& acc) const { acc.addPoint(this->currentPoint(),currentWeight()); }

// ACCESS
  /** Retrieve the current weight
   * \return The value of the weight of the current point
   * \note Valid only if hasWeight() is true
   * \sa hasWeight()
   */
  Scalar currentWeight() const { return mrWeight(1,BaseType::currentIndex()); }

  /** Retrieve an arbitrary weight
   * \return The value of the point's weight
   * \note Valid only if hasWeight() is true
   * \sa hasWeight()
   */
  Scalar weight(int index) const { return mrWeight(1,index); }

private:
  const DerivedWeight&    mrWeight;
};


// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerPreorderedC_h
