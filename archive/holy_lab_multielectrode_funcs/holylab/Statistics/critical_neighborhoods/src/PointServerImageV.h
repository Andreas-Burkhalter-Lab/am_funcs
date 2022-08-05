#ifndef PointServerImageV_h
#define PointServerImageV_h

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

/**  \brief A PointServer for images that returns only value coordinates (pixel values), served in order of increasing spatial distance from a base grid position.
 *
 * \details #include "PointServerImageV.h" <BR>
 *
 * An image takes one or more values at each pixel.  For example, a grayscale image will have a single value (the intensity), whereas an RGB image will have 3 values for each pixel (the intensity of the red, green, and blue color channels).  Consequently, in general "value" for an image is a vector.\n
 * Points are served in order of increasing spatial distance from a base position, which is "on-grid".  Note we distinguish base position (a particular pixel in an image) from base point (a vector of values).
 * \tparam CoordP The vector type used to represent the spatial position in physical units. This should be a native Eigen type with allocated storage (e.g., VectorXf or Vector3d) and should be floating-point.
 * \tparam CoordV The value type. Even if each pixel is a scalar (e.g., uint8_t), you should specify the type as an Eigen type (e.g., Eigen::Matrix<uint8_t,1,1>).
 */
template <typename CoordP,typename CoordV>
class PointServerImageV
{
public:
  /** The scalar data type for position. Pixel values are cast to this scalar type before being returned. */
  typedef typename CoordP::Scalar Scalar;

  /** The scalar data type for image value, as provided by the image data. */
  typedef typename CoordV::Scalar ScalarV;

  /** A type for storing value vectors cast to CoordP::Scalar type */
  typedef typename Eigen::Matrix<Scalar,CoordV::RowsAtCompileTime,1> CoordVP;

  /** The point type (with allocated storage).  This typename is referenced externally */
  typedef CoordVP PointType;

  /** A type representing spatial position in grid coordinates (integers) */
  typedef typename PointServerGrid<CoordP>::CoordG CoordG;

  /** Return type for value data */
  //typedef typename Eigen::CwiseUnaryOp<Eigen::internal::scalar_cast_op<ScalarV, Scalar>, const Eigen::Map<const CoordV> > CoordVReturn;
  typedef typename Eigen::CwiseUnaryOp<Eigen::internal::scalar_cast_op<ScalarV, Scalar>, const CoordV> CoordVReturn;
  //typedef Eigen::Map<const CoordV> CoordVReturn;
  //typedef CoordVP CoordVReturn;

  /** An expression type for returning arbitrary grid points in physical units*/
  typedef typename PointServerGrid<CoordP>::G2PType G2PType;


// LIFECYCLE
  /** The constructor.
   * \param gn An instance of GridNeighbors
   * \param pImage A pointer to the raw image data
   * \param n_values The number of values per pixel
   * \sa GridNeighbors, PointServerGrid
   * \note You should use basePoint(bp) to set or change the base point
   */
  PointServerImageV(const GridNeighbors<CoordP>& gn,const ScalarV *pImage,int n_values) : mPSG(gn), mpImage(pImage), mNValues(n_values) { init(); restart(); }
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
   * \sa basePoint(), basePositionG()
   */
  void restart() { mPSG.restart(); }

  /** Set the base point (not a spatial position---it's the image value(s))
   * \param bp A reference to a vector of image value(s)
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \note This does not change the base position
   * \sa void basePositionG
   */
  template <typename TV> void basePoint(const Eigen::MatrixBase<TV>& bp) {
    mBasePointV = bp;
    restart();
  }

  /** Set the base position to a particular pixel, in grid coordinates, and use the pixel value at that point as the new base point
   * \param bp A reference to the coordinates of the new base point, in grid units
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint, basePositionIndex
   */
  template <typename T> void basePositionG(const Eigen::MatrixBase<T>& bp) {
    mPSG.basePointG(bp);
    mBasePointV = getValue(mPSG.basePointIndex()).template cast<Scalar>();
    restart();
  }

  /** Set the base position to a particular pixel, using a grid index, and use the pixel value at that pixel as the new base point
   * \param Index The integer index specifying a particular pixel (column-major)
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint, void basePositionG
   */
  void basePositionIndex(int Index) {
    mPSG.basePointIndex(Index);
    mBasePointV = getValue(Index).template cast<Scalar>();
    restart();
  }

  /** Update the basePoint to the mean of a neighborhood
   * \param acc An accumulator
   */
  template <typename Accumulator>
  void update(const Accumulator& acc) {
    basePoint(acc.mean());
  }


// ACCESS
  /** The number of pixels in the image */
  int N() const { return mPSG.N(); }

  /** The number of values associated with each pixel */
  int d() const { return mNValues; }

  /** Return the grid size
   * \return A constant reference to a vector holding the grid size
   */
  const CoordG& sizeSpatial() const { return mPSG.size(); }

  /** Obtain the current position's index.
   * \return The integer index of the current position (zero-offset)
   * \sa currentPoint()
   */
  int currentIndex() const { return mPSG.currentIndex(); }

  /** Obtain the current position, in physical units
   * \return Physical coordinates of the current position
   * \sa currentPoint(), currentPositionG(), currentIndex()
   */
  const typename PointServerGrid<CoordP>::ColPSumType currentPosition() const { return mPSG.currentPoint(); }

  /** Obtain the current position, in grid units
   * \return Grid coordinates of the current position
   * \sa currentPoint(), currentPosition(), currentIndex()
   */
  const typename PointServerGrid<CoordP>::ColGSumType currentPositionG() const { return mPSG.currentPointG(); }

  /** Retrieve the current point
   * \return The image value(s) at the current pixel. Note this is cast to type CoordP::Scalar, to be compatible with an Accumulator.
   * \sa currentIndex()
   */
  CoordVReturn currentPoint() const {
    int cIndex = currentIndex();
    //return getValue(cIndex).template cast<Scalar>();  // FIXME segfault
    mTmpV = getValue(cIndex);
    return mTmpV.template cast<Scalar>();
  }

  /** Retrieve pixel values at an arbitrary pixel
   * \param index The pixels index (column-major ordering)
   * \return The value(s) at the indicated pixel
   * \sa currentPoint()
   */
  CoordVReturn point(int index) const {
    //return getValue(index).template cast<Scalar>();  // FIXME segfault
    mTmpV = getValue(index);
    return mTmpV.template cast<Scalar>();
  }

  /** Retrieve the base point
   * \return The value(s) of the base point
   */
  const CoordVP& basePoint() const {
    return mBasePointV;
  }

  /** Retrieve the grid coordinates of the base position
   * \return A constant reference to the coordinates of the base position, in grid units
   */
  const CoordG& basePositionG() const { return mPSG.basePointG(); }

  /** Retrieve the grid index of the base position
   * \return The integer index of the on-grid base position (column-major ordering)
   */
  int basePositionIndex() const { return mPSG.basePointIndex(); }


// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence.
   * \return True (boolean) if at the end.
   */
  bool isAtEnd() const { return mPSG.isAtEnd(); }

  /** Test whether the current point "ties" the previous one, meaning it has the same distance.
   * \return True (boolean) if it is a tie.
   */
  bool isTiesPrevious() const { return mPSG.isTiesPrevious(); }


protected:
private:
  typedef typename Eigen::RowVectorXi RowVectori;
  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
  enum { NeedsToAlign = (sizeof(CoordVP)%16) == 0};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  PointServerGrid<CoordP> mPSG;  // the grid PointServer
  const ScalarV*          mpImage;   // pointer to the image data
  RowVectori mValueOffset;  // relative memory offset for start of nbr values
  CoordVP mBasePointV;      // base point, value coordinates (cast->Scalar)
  int mNValues;   // # of values per pixel
  mutable CoordVP mTmpVP;
  mutable CoordV mTmpV;

  void init() {
    // Store the neighbor offsets for pixel values
    mValueOffset = mPSG.neighborOffsetIndex() * mNValues;
    // Look up the corresponding value for the default base position
    mBasePointV = getValue(mPSG.basePointIndex()).template cast<Scalar>();
  }

  typename Eigen::Map<const CoordV> getValue(int index) const {
    return Eigen::Map<const CoordV>(mpImage + index*mNValues, mNValues);
  }

  void nextPoint() { mPSG++; }


  // Prevent copying
  PointServerImageV(const PointServerImageV& from);
  PointServerImageV& operator=(const PointServerImageV& from);
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerImageV_h

