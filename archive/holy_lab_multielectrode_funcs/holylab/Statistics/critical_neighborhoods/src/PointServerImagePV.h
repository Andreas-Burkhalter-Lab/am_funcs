#ifndef PointServerImagePV_h
#define PointServerImagePV_h

// SYSTEM INCLUDES
//
#include <stdint.h>

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "PointServerGrid.h"
#include "Metric.h"

// FORWARD REFERENCES
//

/** \class ImagePointRef
 * \brief A type that representions location and value for pixels in an image. 
 * \details This class exists so that PointServerImagePV can return a single "point" that really has two components, the spatial location and the value(s) measured at a particular pixel. Typically this "wrapper" will be optimized away by the compiler.\n
 * This version is used to return ImagePoints with class-internal storage, e.g., the basepoint.  If you want to return an expression, see ImagePoint.
 * \tparam PositionType The vector type used to represent the spatial position.
 * \tparam ValueType The value type. Even if it's a scalar (e.g., grayscale images), this will typically be an Eigen type (e.g., Matrix<uint16_t, 1, 1>).
 * \sa ImagePoint
 */
template <typename PositionType,typename ValueType>
class ImagePointRef
{
public:
// LIFECYCLE
  /** Constructor
   * \params p A vector (or expression or map or...) specifying spatial position
   * \params v A scalar (for grayscale images) or vector (for, e.g., RGB images) specifying value
   */
  ImagePointRef(PositionType &p,ValueType &v) : mrPosition(p), mrValue(v) {;}
  ImagePointRef(const ImagePointRef& ip) : mrPosition(ip.position()), mrValue(ip.value()) {;}
  // Default destructor
// OPERATORS
  ImagePointRef& operator=(const ImagePointRef &ip) {
    mrPosition = ip.position();
    mrValue = ip.value();
  }

// ACCESS
  /** Extract the position component
   * \return The position */
  const PositionType& position() const { return mrPosition; }

  /** Extract the value component
   * \return The value */
  const ValueType& value() const { return mrValue; }


private:
  PositionType& mrPosition;
  ValueType&    mrValue;
};


/** \class ImagePoint
 * \brief A type that representions location and value for pixels in an image. 
 * \details This class exists so that PointServerImagePV can return a single "point" that really has two components, the spatial location and the value(s) measured at a particular pixel. Typically this "wrapper" will be optimized away by the compiler.\n
 * This class is the one to use if you want to return an expression.  If instead you want to return a position/value pair with class-internal storage (e.g., the basepoint), use ImagePointRef.
 * \tparam PositionType The vector type used to represent the spatial position.
 * \tparam ValueType The value type. Even if it's a scalar (e.g., grayscale images), this will typically be an Eigen type (e.g., Matrix<uint16_t, 1, 1>).
 * \sa ImagePointRef
 */
template <typename PositionType,typename ValueType>
class ImagePoint
{
public:
// LIFECYCLE
  /** Constructor
   * \params p A vector (or expression or map or...) specifying spatial position
   * \params v A scalar (for grayscale images) or vector (for, e.g., RGB images) specifying value
   */
  ImagePoint(const PositionType &p,const ValueType &v) : mrPosition(p), mrValue(v) {;}
  template <typename PT,typename VT>
  ImagePoint(const ImagePoint<PT,VT>& ip) : mrPosition(ip.position()), mrValue(ip.value()) {;}
  // Default destructor
// OPERATORS
  template <typename PT,typename VT>
  ImagePoint& operator=(const ImagePoint<PT,VT>& ip) {
    mrPosition = ip.position();
    mrValue = ip.value();
  }
// ACCESS
  /** Extract the position component
   * \return The position */
  const PositionType& position() const { return mrPosition; }

  /** Extract the value component
   * \return The value */
  const ValueType& value() const { return mrValue; }


private:
  PositionType mrPosition;
  ValueType    mrValue;
};


/**  \brief A PointServer for images, mixing spatial and value coordinates to define total distance
 *
 * \details #include "PointServerImagePV.h" <BR>
 *
 * An image takes one or more values at each pixel.  For example, a grayscale image will have a single value (the intensity), whereas an RGB image will have 3 values for each pixel (the intensity of the red, green, and blue color channels).  Consequently, in general "value" for an image is a vector.\n
 * Points are served in order of increasing total distance.
 * The total distance \f$d\f$ between the base point \f$({\bf x}_0,{\bf v}_0)\f$ (position \f${\bf x}_0\f$, with pixel value(s) \f${\bf v}_0\f$) and a different pixel \f$({\bf x},{\bf v})\f$ is defined as
 * \f[ d^2 = |{\bf x} - {\bf x}_0|^2 + c^2 d_v^2({\bf v}_0,{\bf v}) \f]
 * where \f$d_v\f$ is the "value metric" distance between the pixel values, and \f$c\f$ is a coefficient which expresses the tradeoff between spatial and value distances.
 *
 * \tparam CoordP The vector type used to represent the spatial position, in physical units. This should be a native Eigen type with allocated storage (e.g., VectorXf or Vector3d) and should probably be floating-point.
 * \tparam CoordV The value type. Even if each pixel is a scalar (e.g., uint8_t), you should specify the type as an Eigen type (e.g., Eigen::Matrix<uint8_t,1,1>).
 * \tparam MetricV The metric (see Metric.h) used to measure image value distances.  Default is Euclidean. Note that while the metric operates on value coordinates, to prevent overflow (e.g., for uint8 images) the value coordinates will first be converted to the same scalar type as used in CoordP.  Make sure you declare your metric appropriately.
 * \tparam updateValueCoefficient If true, when updating the value coefficient is set to balance the positional and value variance from the previous neighborhood.  If false, the value coefficient is held constant.
 *  
 */
template <typename CoordP,typename CoordV,typename MetricType = typename Metric::Euclidean,bool updateValueCoefficient = false>
class PointServerImagePV
{
public:
  /** The scalar data type for position. This is the type in which combined distance will be measured. */
  typedef typename CoordP::Scalar Scalar;

  /** The scalar data type for image value. */
  typedef typename CoordV::Scalar ScalarV;

  /** A type for storing value vectors cast to CoordP::Scalar type */
  typedef typename Eigen::Matrix<Scalar,CoordV::RowsAtCompileTime,1> CoordVP;

  /** A type representing spatial position in grid coordinates (integers) */
  typedef typename PointServerGrid<CoordP>::CoordG CoordG;

  /** Return type for position data */
  typedef typename PointServerGrid<CoordP>::ColPSumType CoordPReturn;

  /** Return type for value data */
  typedef typename Eigen::Map<const CoordV> CoordVReturn;

  /** An expression type for returning arbitrary grid points in physical units*/
  typedef typename PointServerGrid<CoordP>::G2PType G2PType;


// LIFECYCLE
  /** The constructor.
   * \param gn An instance of GridNeighbors
   * \param pImage A pointer to the raw image data
   * \param n_values The number of values per pixel
   * \param metric An instance of a metric class
   * \sa GridNeighbors, PointServerGrid
   * \note You should use basePoint(bp) to set or change the base point
   * \note The default ValueCoefficient is zero, you probably want to set this first
   */
  PointServerImagePV(const GridNeighbors<CoordP>& gn,const ScalarV *pImage,int n_values, const MetricType &metric = Metric::Euclidean()) : mPSG(gn), mpImage(pImage), mNValues(n_values), mSDistance(mPSG.nNeighbors()), comp(mSDistance), mrMetric(metric) { init(); restart(); }
  // Prevent copying: private copy constructor & assignment operator
  // Default destructor


// OPERATORS
  /** Prefix increment operator
   * \return Nothing
   */
  void operator++() { nextPoint(); }

  /** Postfix increment operator
   * \return Nothing
   */
  void operator++(int unused) { nextPoint(); }
  

// OPERATIONS
  /** Restart at the beginning, keeping the current base point.
   * \sa basePoint(), basePointG()
   */
  void restart() {
    mPSG.restart();
    mSDistanceLast = -1;
    mGridNbrIndex.clear();
    nextPoint();
  }

  /** Set the base point to a position/value pair
   * \param bp A reference to an ImagePoint
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePointG
   */
  template <typename TP,typename TV> void basePoint(const ImagePoint<TP,TV>& bp) {
    mPSG.basePoint(bp.position());  // this stores a grid-rounded basepoint
    mBasePointP = bp.position();    // store the actual position basepoint
    mBasePointV = bp.value();
    restart();
  }

  // Also a version for ImagePointRef
  template <typename TP,typename TV> void basePoint(const ImagePointRef<TP,TV>& bp) {
    mPSG.basePoint(bp.position());  // this stores a grid-rounded basepoint
    mBasePointP = bp.position();    // store the actual position basepoint
    mBasePointV = bp.value();
    restart();
  }

  /** Set the base point to an on-grid point
   * This is used when "starting fresh" at a particular pixel
   * \param bp A reference to the coordinates of the new base point, in grid units
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint
   */
  template <typename T> void basePointG(const Eigen::MatrixBase<T>& bp) {
    mPSG.basePointG(bp);
    mBasePointP = mPSG.basePoint();
    mBasePointV = getValue(mPSG.basePointIndex()).template cast<Scalar>();
    restart();
  }

  /** Set the base point to an on-grid point, using a grid index
   * This is used when "starting fresh" at a particular pixel
   * \param Index The integer index specifying a particular grid point (column-major)
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint, void basePointG
   */
  void basePointIndex(int Index) {
    mPSG.basePointIndex(Index);
    mBasePointP = mPSG.basePoint();
    mBasePointV = getValue(mPSG.basePointIndex()).template cast<Scalar>();
    restart();
  }

  /** Set the coefficient setting the balance between value and position
   * "Units" are "position/value"
   * \param v2p The coefficient converting value into position
   */
  void valueCoefficient(Scalar v2p) { mCoeffValue = v2p*v2p; }

  /** Set the value coefficient based on the ratio of variances in the first k "shells."
   * \param k The number of "shells" over which to calculate variance
   * \note This calls restart().
   * \note If there is no variance in value over the first k shells, more shells will be added until some non-zero variance is obtained (or all neighboring points exhausted).
   */
  void balanceValueCoefficient(int k) {
    int nNonTies;
    Scalar cp,cv;

    nNonTies = 0;
    cp = cv = 0;
    mPSG.restart();
    while ((nNonTies < k || cv == 0) && !mPSG.isAtEnd()) {
      mTmpV = getValue(mPSG.currentIndex()).template cast<Scalar>() - mBasePointV;
      cv += mTmpV.squaredNorm();
      cp += mPSG.currentOffset().squaredNorm();
      mPSG++;
      nNonTies += !mPSG.isTiesPrevious();
    }
    if (cv > 0)
      mCoeffValue = cp/cv;
    restart();
  }

  /** Update the basePoint and (optionally) the value coefficient
   * \param acc An accumulator of type AccumulatorImagePV
   *
   */
  template <typename AccumulatorType>
  void update(const AccumulatorType& acc) {
    basePoint(acc.mean());

    if (updateValueCoefficient) {
      Scalar cv = acc.value().covariance_n();
      Scalar cp = acc.position().covariance_n();
      Scalar coeffValueNew;
      if (cv > 0)
	coeffValueNew = cp/cv;
      else
	coeffValueNew = mCoeffValue;
      // Limit the change to twofold on each iteration
      if (coeffValueNew < mCoeffValue/2)
	coeffValueNew = mCoeffValue/2;
      if (coeffValueNew > 2*mCoeffValue)
	coeffValueNew = 2*mCoeffValue;
      mCoeffValue = coeffValueNew;
    }
  }
  

// ACCESS
  /** Return the number of pixels in the image
   * \return An integer with the number of pixels in the image
   */
  int N() const { return mPSG.N(); }

  /** Return the grid size
   * \return A constant reference to a vector holding the grid size
   */
  const CoordG& sizeSpatial() const { return mPSG.size(); }

  /** Obtain the current point's index.
   * \return The integer index of the current point (zero-offset)
   * \sa currentPoint()
   */
  int currentIndex() const { return mGridIndex[mGridNbrIndex[0]]; }

  /** Retrieve the current point
   * \return An ImagePoint containing the spatial (physical units) and value coordinates of the current point
   * \sa currentIndex()
   */
  ImagePoint<const CoordPReturn,const CoordVReturn> currentPoint() const {
    int cIndex = currentIndex();
    return ImagePoint<const CoordPReturn,const CoordVReturn>(mPSG.neighborPoint(mGridNbrIndex[0]),getValue(cIndex));
  }

  /** Retrieve the grid coordinates of the current point
   * \return A reference to the grid coordinates
   * \sa currentIndex()
   */
  const typename PointServerGrid<CoordP>::ColGSumType currentPointG() const {
    return mPSG.neighborPointG(mGridNbrIndex[0]);
  }

  /** Retrieve an arbitrary on-grid point
   * \param index The on-grid point's index (column-major)
   * \return An ImagePoint containing the spatial (in physical units) and value coordinates of the point.
   * \sa currentPoint()
   */
  //ImagePoint<const G2PType,const CoordVReturn> point(int index) const {
  //  return ImagePoint<const G2PType,const CoordVReturn>(mPSG.point(index),getValue(index));
  ImagePointRef<const CoordP,const CoordV> point(int index) const {
    mTmpV = getValue(index);
    return ImagePointRef<const CoordP,const CoordV>(mPSG.point(index),mTmpV);
  }

  /** Retrieve the base point
   * \return An ImagePoint containing the spatial and value coordinates of the base point
   */
  ImagePointRef<const CoordP,const CoordVP> basePoint() const {
    return ImagePointRef<const CoordP,const CoordVP>(mBasePointP,mBasePointV);
  }

  /** Retrieve the grid coordinates of the base point
   * \return A constant reference to the coordinates of the base point, in grid units
   */
  const CoordG& basePointG() const { return mPSG.basePointG(); }

  /** Retrieve the grid index of the base point
   * \return The integer index of the on-grid base point (column-major ordering)
   */
  int basePointIndex() const { return mPSG.basePointIndex(); }

  /** Get the coefficient setting the balance between value and position
   * \return The coefficient converting value to position
   */
  Scalar valueCoefficient() const { return std::sqrt(mCoeffValue); }


// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence.
   * \return True (boolean) if at the end.
   */
  bool isAtEnd() const { return (mGridNbrIndex.empty() && mPSG.isAtEnd()); }

  /** Test whether the current point "ties" the previous one, meaning it has the same distance.
   * \return True (boolean) if it is a tie.
   */
  bool isTiesPrevious() const { 
    if (mGridNbrIndex.empty())
      return false;
    else
      return (mSDistance[mGridNbrIndex[0]] == mSDistanceLast);
  }

  /** Report the size of the neighbor heap (useful for diagnostics)
   * \return The number of points currently on the heap
   */
  int heapSize() const { return mGridNbrIndex.size(); }

protected:
private:
  typedef typename Eigen::Matrix<Scalar,1,Eigen::Dynamic> RowVectorP;
  typedef typename Eigen::RowVectorXi RowVectori;
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
  enum { NeedsToAlign = ((sizeof(CoordP)%16)==0 ||
			 (sizeof(CoordV)%16)==0 ||
			 (sizeof(CoordVP)%16) == 0)};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  PointServerGrid<CoordP> mPSG;  // the grid PointServer
  const ScalarV*          mpImage;   // pointer to the image data
  RowVectorP  mSDistance;    // full distance to neighbors, indexed by nbrIndex
  IndirectGreater<RowVectorP> comp;  // comparison class for heap
  const MetricType&       mrMetric;    // reference to metric object
  std::vector<int>        mGridNbrIndex;  // the heap of grid neighbor indices
  RowVectori mGridIndex;     // grid index of neighbors, indexed by nbrIndex
  RowVectori mValueOffset;  // relative memory offset for start of nbr values
  CoordP mBasePointP;       // base point, spatial coordinates (phys. units)
  CoordVP mBasePointV;      // base point, value coordinates (cast->Scalar)
  Scalar mCoeffValue;       // the coefficient converting value^2->pos^2
  Scalar mSDistanceLast;    // full distance to previous neighbor
  mutable CoordVP mTmpV;  // temporary storage of image values cast to Scalar
  int mNValues;   // # of values per pixel

  void init() {
    // Allocate storage
    int d = mPSG.size().size();
    int n = mPSG.nNeighbors();
    int i;
    //std::cout << "d = " << d << ", n = " << n << std::endl;
    mGridNbrIndex.reserve(n);
    mGridIndex.resize(n);
    mValueOffset = mPSG.neighborOffsetIndex() * mNValues;
    //std::cout << "mValueOffset: " << mValueOffset << std::endl;
    mBasePointP.resize(d,1);
    mTmpV.resize(mNValues,1);
    // Look up the corresponding value for the default base point
    mBasePointV = getValue(mPSG.basePointIndex()).template cast<Scalar>();
    //std::cout << "mBasePointV " << mBasePointV.transpose() << std::endl;
    mCoeffValue = 0;
  }

  CoordVReturn getValue(int index) const {
    return CoordVReturn(mpImage + index*mNValues, mNValues);
  }

  void nextPoint() {
    Scalar currentSDistance,thisSDistance;
    int thisNbrIndex, thisIndex;
    if (!mGridNbrIndex.empty()) {
      mSDistanceLast = mSDistance[mGridNbrIndex[0]];
      pop_heap(mGridNbrIndex.begin(),mGridNbrIndex.end(),comp);
      mGridNbrIndex.erase(mGridNbrIndex.end()-1);
    }
    if (!mGridNbrIndex.empty())
      currentSDistance = mSDistance[mGridNbrIndex[0]];
    else
      currentSDistance = std::numeric_limits<Scalar>::max();  // sentinel
    /*
    std::cout << "PSIPV::nextPoint:\ncurrentSDistance " << currentSDistance
	      << ", mPSG.curSD() " << mPSG.currentSquareDist()
	      << ", mPSG.isAtEnd() " << mPSG.isAtEnd() << std::endl;
    */
    while (!mPSG.isAtEnd() && mPSG.currentSquareDist() <= currentSDistance) {
      //std::cout << "curSD = " << currentSDistance << ", mPSG.curSD = " << mPSG.currentSquareDist() << std::endl;
      // Can't be sure that we know the closest neighbor, must add more points
      thisNbrIndex = mPSG.currentNeighborIndex();
      thisIndex = mPSG.currentIndex();
      //std::cout << "mPSG: currentPoint " << mPSG.currentPoint().transpose() << ", thisNbrIndex " << thisNbrIndex << ", thisIndex " << thisIndex << std::endl;
      // Value component
      mTmpV = getValue(thisIndex).template cast<Scalar>();
      //std::cout << "mTmpV " << mTmpV.transpose() << std::endl;
      thisSDistance = mrMetric.squaredDistance(mBasePointV,mTmpV);
      thisSDistance *= mCoeffValue;
      // Spatial component
      thisSDistance += (mBasePointP - mPSG.currentPoint()).squaredNorm();
      if (!std::isfinite(thisSDistance))
	throw std::runtime_error("Combined distance is not finite");
      //std::cout << "Total sdist " << thisSDistance << std::endl;
      // Store the relevant data for this neighbor
      mSDistance[thisNbrIndex] = thisSDistance;
      mGridIndex[thisNbrIndex] = thisIndex;
      // Add to heap
      mGridNbrIndex.push_back(thisNbrIndex);
      push_heap(mGridNbrIndex.begin(),mGridNbrIndex.end(),comp);
      //std::cout << "Heap: ";
      //for (int i = 0; i < mGridNbrIndex.size(); i++)
      //	std::cout << mGridNbrIndex[i] << ' ';
      //std::cout << std::endl;
      // Advance to the next grid point
      mPSG++;
      // For the next comparison, use the closest point so far
      currentSDistance = mSDistance[mGridNbrIndex[0]];
    }
    //std::cout << "ALL DONE WITH nextPoint()" << std::endl;
  }

  // Prevent copying
  PointServerImagePV(const PointServerImagePV& from);
  PointServerImagePV& operator=(const PointServerImagePV& from);
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerImagePV_h

