#ifndef PointServerGrid_h
#define PointServerGrid_h

/** \class PointServerGrid
 *  \brief A PointServer ordering points by distance from a base point on a grid
 *
 * \details #include "PointServerGrid.h" <BR>
 *
 * This serves points in order of increasing distance from a base point.  Since we can use the grid geometry to identify neighbors, this is much more efficient than checking each point.\n
 * The main subtlety with this class is that there are two notions of location.   "Grid coordinates" are integers that identify particular grid points.  "Physical coordinates" are positions in space, measured in physical-space units. In the simplest case physical space is a scaled version of grid space, with each coordinate independently scaled. In more complex cases, physical space can be a skewed version of the grid.\n
 * To avoid the need for data duplication when there are multiple PointServers for a single grid geometry (e.g., when multithreading), the storage for all the distance-ranked grid displacements is split out in a separate GridNeighbors structure.
 * \tparam CoordP The vector type used to represent the spatial position, in physical units. This should be a native Eigen type with allocated storage (e.g., Vector2d or VectorXf).
 * \sa GridNeighbors
 *  
 */

// SYSTEM INCLUDES
//
#include <functional>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "Eigen/Dense"
#include "Eigen/SVD"

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "IndirectCompare.h"

// FORWARD REFERENCES
//

/** \class GridNeighbors
 * \brief A structure for storing displacements to neighboring grid points in order of increasing physical distance.
 * \details Multiple PointServerGrid instances can be initialized using a single instance of GridNeighbors, avoiding the need to duplicate this information (which can be sizable, if the maximum radius is large).  To make access simple, the data members are public.
 * \tparam CoordP The vector type used to represent the spatial position, in physical units. This should be a native Eigen type with allocated storage (e.g., Vector2d or VectorXf).
 * \sa PointServerGrid
 */
template <typename CoordP>
class GridNeighbors
{
public:
  /** The scalar data type */
  typedef typename CoordP::Scalar Scalar;

  /** A type representing spatial position in grid coordinates (integers) */
  typedef typename Eigen::Matrix<int,CoordP::RowsAtCompileTime,1> CoordG;

  /** An "array" of points/displacements, in physical units */
  typedef typename Eigen::Matrix<Scalar,CoordP::RowsAtCompileTime,Eigen::Dynamic> MatrixP;

  /** An "array" of points/displacements, in grid units */
  typedef typename Eigen::Matrix<int,CoordP::RowsAtCompileTime,Eigen::Dynamic> MatrixG;

  /** The type storing square distances of grid displacements */
  typedef typename Eigen::Matrix<Scalar,1,Eigen::Dynamic> RowVectorP;

  /** The transform type for converting grid coords to physical coords */
  typedef typename Eigen::Matrix<Scalar,CoordP::RowsAtCompileTime,CoordP::RowsAtCompileTime> TformType;

  typedef typename Eigen::JacobiSVD<TformType> DecompType;

  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  enum { NeedsToAlign = ((sizeof(CoordP)%16)==0 ||
			 (sizeof(CoordG)%16)==0)};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
#endif // DOXYGEN_SHOULD_SKIP_THIS

  CoordG     mSz;        // the grid size
  CoordG     mCumSz;     // per-coordinate grid index (column-major)
  TformType  mG2P;       // the matrix mapping grid coords to phys coords
  DecompType mG2Pdecomp; // the decomposition of mG2P (for inversion)
  TformType  mG2Pinv;    // the inverse mapping
  Scalar     mR;         // the size of the neighborhood radius (phys coords)
  MatrixP    mOffsetP;   // offsets of all neighboring grid points (phys coords)
  MatrixG    mOffsetG;   // offsets of all neighboring grid points (grid coords)
  Eigen::RowVectorXi   mOffsetIndex;   // index offset of neighboring grid pts
  RowVectorP mSDistance; // square distance to neighboring grid points (phys)

// LIFECYCLE
  /** The constructor.
   * \param sz A vector specifying the size of the grid in each dimension (i.e., the number of grid points along each coordinate).
   * \param r The radius of the local neighborhood, in physical units
   * \param g2p An optional square matrix mapping grid coordinates (g) to physical coordinates (p). The default is the identity, in which case grid coordinates are the same as physical coordinates.
   */
  template <typename T1,typename T2>
  GridNeighbors(const Eigen::MatrixBase<T1> &sz,Scalar r,const Eigen::MatrixBase<T2> &g2p) : mSz(sz), mG2P(g2p), mR(r) { init(); }
  template <typename T1>
  GridNeighbors(const Eigen::MatrixBase<T1> &sz,Scalar r) : mSz(sz), mG2P(MatrixP::Identity(sz.size(),sz.size())), mR(r) { init(); }
  // Prevent copying: private copy constructor & assignment operator
  // Default destructor


// ACCESS
  /** Return the number of points in the grid
   * \return An integer with the number of points in the grid
   */
  int N() const { return mCumSz[mCumSz.size()-1]*mSz[mSz.size()-1]; }


// INQUIRY
  /** \return The maximum number of "tying" points */
  int maxTyingPoints() const {
    int i,n,nmax;
    n = nmax = 0;
    for (i = 1; i < mSDistance.size(); i++) {
      if (mSDistance(i) == mSDistance(i-1))
	n++;
      else {
	if (n > nmax)
	  nmax = n;
	n = 0;
      }
    }
    return nmax+1;
  }

private:
  // Prevent copying
  GridNeighbors(const GridNeighbors& from);
  GridNeighbors& operator=(const GridNeighbors& from);

  void init() {
    int i;
    // Allocate storage
    int d = mSz.size();
    if (mG2P.rows() != d || mG2P.cols() != d)
      throw std::runtime_error("Dimensionality is inconsistent");
    mCumSz.resize(d,1);
    // Determine grid size indexing
    mCumSz[0] = 1;
    for (i = 1; i < d; i++)
      mCumSz[i] = mCumSz[i-1]*mSz[i-1];
    // Compute the decomposition of the transform matrix
    mG2Pdecomp.compute(mG2P,Eigen::ComputeFullU | Eigen::ComputeFullV);
    if ((mG2Pdecomp.singularValues().array() == 0).any())
      throw std::runtime_error("The transform matrix is degenerate");
    // Compute the inverse of the decomposition (calculate ahead
    // because we will need it many times, and this avoids calls to
    // posix_memalign)
    mG2Pinv = mG2Pdecomp.matrixV()
      * mG2Pdecomp.singularValues().cwiseInverse().asDiagonal()
      * mG2Pdecomp.matrixU().adjoint();
    // Calculate the maximum displacement along each grid coordinate
    // that fits within a sphere of radius mR
    // (Use the SVD: the column vectors of V, divided by the
    // corresponding singular value, give the extent along each
    // coordinate)
    CoordP maxD = ((mG2Pdecomp.matrixV()*
		   mG2Pdecomp.singularValues().cwiseInverse().asDiagonal())
		   .array().abs().rowwise().maxCoeff()) * mR;
    CoordG maxDG = maxD.unaryExpr(std::ptr_fun<Scalar,Scalar>(std::ceil)).template cast<int>();
    maxDG = maxDG.cwiseMin(mSz); // no displacement bigger than the grid size needs to be considered
    // Generate all possible displacements
    int n = (2*maxDG.array()+1).prod();   // number of vector displacements to consider
    MatrixP testOffsetG(d,n);     // Scalar type so can compute matrix product
    testOffsetG.col(0) = -maxDG.template cast<Scalar>();  // initialize at left edge
    int coordindx;
    Scalar *pOffset;
    for (pOffset = testOffsetG.data()+d; pOffset < testOffsetG.data()+d*n; pOffset+=d) {
      coordindx = 0;
      // increment-with-carry
      while ((pOffset[coordindx] = pOffset[coordindx-d]+1) > maxDG[coordindx]) {
	pOffset[coordindx] = -maxDG[coordindx];
	coordindx++;
	if (coordindx >= d)
	  break;
      }
      //assert(coordindx <= d || pOffset+d >= testOffsetG.data()+d*n);
      // copy forward the ones not involved in the carry
      for (coordindx++; coordindx < d; coordindx++)
	pOffset[coordindx] = pOffset[coordindx-d];
    }
    MatrixP testOffsetP = mG2P*testOffsetG;
    RowVectorP testSDistance = testOffsetP.colwise().squaredNorm();
    // To avoid being confused by roundoff error, standardize the accuracy of each square distance
    Scalar mR2 = mR*mR;
    Scalar resolution = 10*std::numeric_limits<Scalar>::epsilon()*mR2;
    Eigen::RowVectorXi testSDistanceI = (testSDistance.array()/resolution + 0.5).template cast<int>();
    testSDistance = resolution*testSDistanceI.template cast<Scalar>();
    // Sort them in order of increasing square distance
    std::vector<int> sortOrder(n);
    for (i = 0; i < n; i++)
      sortOrder[i] = i;
    IndirectLess<RowVectorP> comp(testSDistance);
    std::sort(sortOrder.begin(),sortOrder.end(),comp);
    // Find the cutoff for where the distance becomes larger than mR
    for (i = 0; i < n; i++)
      if (testSDistance[sortOrder[i]] > mR2)
	break;
    n = i;
    // Keep those points that are within radius mR
    mOffsetP.resize(d,n);
    mOffsetG.resize(d,n);
    mSDistance.resize(1,n);
    mOffsetIndex.resize(n);
    for (i = 0; i < n; i++) {
      int j = sortOrder[i];
      mOffsetP.col(i) = testOffsetP.col(j);
      mOffsetG.col(i) = testOffsetG.col(j).template cast<int>();
      mOffsetIndex(i) = (mOffsetG.col(i).array() * mCumSz.array()).sum();
      mSDistance(i) = testSDistance(j);
    }
  }
};

template <typename CoordP>
class PointServerGrid
{
public:
  /** The scalar data type */
  typedef typename CoordP::Scalar Scalar;

  /** A type representing spatial position in grid coordinates (integers) */
  typedef typename Eigen::Matrix<int,CoordP::RowsAtCompileTime,1> CoordG;

  /** The expression type of currentPoint() (when evaluated, yields a column vector<Scalar>) */
  typedef typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<Scalar>, typename GridNeighbors<CoordP>::MatrixP::ConstColXpr,const CoordP> ColPSumType;

  /** The expression type of currentPointG() (when evaluated, yields a column vector<int>) */
  typedef typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<int>, typename GridNeighbors<CoordP>::MatrixG::ConstColXpr,const CoordG> ColGSumType;

  /** The return expression type of point(index) */
  typedef typename Eigen::ProductReturnType<const typename GridNeighbors<CoordP>::TformType,const CoordP>::Type G2PTypeConst;
  typedef typename Eigen::ProductReturnType<typename GridNeighbors<CoordP>::TformType,CoordP>::Type G2PType;

// LIFECYCLE
  /** The constructor.
   * \param gn An instance of GridNeighbors.
   * \note You should use basePoint(bp) to set or change the base point; the PointServer is not ready to use until you set the base point
   */
  PointServerGrid(const GridNeighbors<CoordP>& gn) : mGN(gn) { init(); restart(); }
  // Default copy constructor & assignment operator
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
  void restart() { mIndex = 0; mSDistanceLast = -1; }

  /** Set the base point, in physical units
   * \param bp A reference to the coordinates of the new base point, in physical units
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \note bp does not have to be on-grid, but it will be "rounded off" to the nearest grid point. Consequently, if you call basePoint() you may not get exactly bp back.
   * \sa void basePointG
   */
  template <typename T> void basePoint(const Eigen::MatrixBase<T>& bp) {
    // Convert physical coordinates to grid coordinates
    //mTmpP2 = mGN.mG2Pdecomp.solve(bp);  // "inverse" of mG2P
    mTmpP = bp;
    mTmpP2.noalias() = mGN.mG2Pinv*mTmpP;
    mBasePointG = (mTmpP2.array()+0.5).unaryExpr(std::ptr_fun<Scalar,Scalar>(std::floor)).template cast<int>();  // round off
    basePointBounds();
    // Set the actual physical coordinates to coincide with the grid coordinates
    mTmpP = mBasePointG.template cast<Scalar>();
    mBasePointP.noalias() = mGN.mG2P*mTmpP;
    restart();
  }

  /** Set the base point, in grid units
   * \param bp A reference to the coordinates of the new base point, in grid units
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint
   */
  template <typename T> void basePointG(const Eigen::MatrixBase<T>& bp) {
    mBasePointG = bp;
    basePointBounds();
    mTmpP = mBasePointG.template cast<Scalar>();
    mBasePointP.noalias() = mGN.mG2P*mTmpP;
    restart();
  }

  /** Set the base point in terms of a grid index
   * \param Index The integer index specifying a particular grid point (column-major)
   * \exception std::runtime_error if you specify a base point that is not within the confines of the grid
   * \sa void basePoint
   */
  void basePointIndex(int Index) {
    index2grid(Index);
    mBasePointG = mTmpP.template cast<int>();
    basePointBounds();
    mBasePointP.noalias() = mGN.mG2P*mBasePointG.template cast<Scalar>();
    restart();
  }    


// ACCESS
  /** Return the number of points in the grid
   * \return An integer with the number of points in the grid
   */
  int N() const { return mGN.N(); }

  /** Return the grid size
   * \return A constant reference to a vector holding the grid size
   */
  const CoordG& size() const { return mGN.mSz; }

  /** Return the number of potential neighbors within the radius specified at construction time.
   * Note that the actual number of neighbors may be smaller if the base point is near the edge of the grid.
   * \return An integer nnbrs, the number of potential neighbors
   */
  int nNeighbors() const { return mGN.mOffsetIndex.size(); }

  /** Return the matrix of all grid offsets within the radius specified at construction time. Each column is a separate offset.
   * \return A ndims-by-nnbrs matrix of integers
   */
  const typename GridNeighbors<CoordP>::MatrixG& neighborOffsetG() const { return mGN.mOffsetG; }

  /** Return the row vector of all grid index offsets within the radius specified at construction time.
   * \return A row vector of length nnbrs
   */
  const Eigen::RowVectorXi& neighborOffsetIndex() const { return mGN.mOffsetIndex; }

  /** Obtain the current point's ranking among neighbors.
   * \return The integer index of the current neighbor (zero-offset)
   */
  int currentNeighborIndex() const { return mIndex; }

  /** Obtain the current point's grid index
   * \return The integer index of the current point (zero-offset) in the grid
   */
  int currentIndex() const { return mBasePointIndex+mGN.mOffsetIndex[mIndex]; }

  /** Retrieve the current grid point in physical coordinates.
   * The current point incorporates the base point's position.
   * \return An expression containing the physical coordinates of the current grid point
   * \sa currentPoint()
   */
  // Note: return as an expression (basepoint can be
  // re-subtracted as a no-op?)
  const ColPSumType currentPoint() const {
    return mGN.mOffsetP.col(mIndex) + mBasePointP;
  }

  /** Retrieve an arbitrary grid point in physical coordinates.
   * \param Index The absolute grid index of the point
   * \return A const reference to the physical coordinates of the grid point
   * \sa currentPoint()
   */
  //  const G2PType point(int Index) const {
  const CoordP& point(int Index) const {
    index2grid(Index);
    // Convert to physical coordinates
    //return mGN.mG2P*mTmpP;
    mTmpP2.noalias() = mGN.mG2P*mTmpP;
    return mTmpP2;
  }

  /** Retrieve a neighbor grid point in physical coordinates.
   * \param nbrIndex The neighbor index of the point
   * \return An expression containing the physical coordinates of the current grid point
   * \sa currentPoint()
   */
  const ColPSumType neighborPoint(int nbrIndex) const {
    return mGN.mOffsetP.col(nbrIndex) + mBasePointP;
  }

  /** Retrieve a neighbor grid point in grid coordinates.
   * \param nbrIndex The neighbor index of the point
   * \return An expression containing the grid coordinates of the current grid point
   * \sa neighborPoint()
   */
  const ColGSumType neighborPointG(int nbrIndex) const {
    return mGN.mOffsetG.col(nbrIndex) + mBasePointG;
  }

  /** Retrieve the current offset in physical coordinates.
   * The offset is relative to the base point.
   * \return A constant reference to the current offset in physical coordinates
   * \sa currentIndex(), currentOffset(), currentPointG()
   */
  typename GridNeighbors<CoordP>::MatrixP::ConstColXpr currentOffset() const {
    return mGN.mOffsetP.col(mIndex);
  }

  /** Retrieve the current grid point in grid coordinates.
   * The current point incorporates the base point's position.
   * \return A constant reference to the grid coordinates of the current grid point
   * \sa currentPoint()
   */
  const ColGSumType currentPointG() const {
    return mGN.mOffsetG.col(mIndex) + mBasePointG;
  }

  /** Retrieve the current offset in grid coordinates.
   * The current offset is relative to the base point.
   * \return A constant reference to the grid coordinates of the current offset
   * \sa currentPoint(), currentOffset()
   */
  typename GridNeighbors<CoordP>::MatrixG::ConstColXpr currentOffsetG() const {
    return mGN.mOffsetG.col(mIndex);
  }

  /** Retrieve the base point, in physical units.
   * \return A constant reference to the coordinates of the base point, in physical units
   * \sa basePointG()
   */
  const CoordP& basePoint() const { return mBasePointP; }

  /** Retrieve the base point, in grid units.
   * \return A constant reference to the coordinates of the base point, in grid units
   * \sa basePoint()
   */
  const CoordG& basePointG() const { return mBasePointG; }

  /** Retrieve the base point's index
   * \return The integer index of the base point (column-major)
   * \sa basePoint()
   */
  int basePointIndex() const { return mBasePointIndex; }

  /** Return the current square-displacement from the base point, in physical units
   */
  Scalar currentSquareDist() { return mGN.mSDistance[mIndex]; }


// INQUIRY
  /** Test whether the PointServer has reached the end of the sequence.  */
  bool isAtEnd() const { return (mIndex >= mGN.mOffsetP.cols()); }

  /** Test whether the current point "ties" the previous one, meaning it has the same distance.  */
  bool isTiesPrevious() const { return (mGN.mSDistance[mIndex] == mSDistanceLast); }

  // The next lines are an Eigen trick to get the compiler to
  // properly align fixed-sized objects
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  enum { NeedsToAlign = ((sizeof(CoordP)%16)==0 ||
			 (sizeof(CoordG)%16)==0)};
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
#endif // DOXYGEN_SHOULD_SKIP_THIS

protected:
private:
  const GridNeighbors<CoordP>& mGN;  // grid displacement "database"
  CoordP     mBasePointP;// current base point, phys coords
  CoordG     mBasePointG;// current base point, grid coords
  int        mBasePointIndex; // current base point index (column-major)
  CoordG     l,u;        // lower and upper bounds on valid displacements
  Scalar     mSDistanceLast;  // storage for evaluating ties
  int        mIndex;     // the current position in the nbr sequence
  mutable CoordP mTmpP,mTmpP2;  // purely temporary storage

  void init() {
    // Allocate storage
    int d = mGN.mSz.size();
    mBasePointP.resize(d,1);
    mBasePointG.resize(d,1);
    l.resize(d,1);
    u.resize(d,1);
    mTmpP.resize(d,1);
    mTmpP2.resize(d,1);
    //  Set default base point
    mBasePointP.setConstant(0);
    mBasePointG.setConstant(0);
    basePointBounds();
  }

  void nextPoint() {
    //typename MatrixG::ColXpr thisColG;
    mSDistanceLast = mGN.mSDistance[mIndex];
    mIndex++;
    const int *pOffsetG;
    int i;
    while (!isAtEnd()) {
      pOffsetG = &mGN.mOffsetG(0,mIndex);  // use pointers as easier/more efficient
      // Is the new point within grid bounds?
      bool ok = true;
      for (i = 0; i < mGN.mOffsetG.rows(); i++, pOffsetG++)
        if ((*pOffsetG < l[i]) || (*pOffsetG >= u[i])) {
	  ok = false;
          break;  
        }
      if (ok)
	break;
      // point is not OK, try the next one
      pOffsetG += mGN.mOffsetG.rows()-i;
      mIndex++;
    }
  }

  void basePointBounds() {
    if ((mBasePointG.array() < 0).any() || (mBasePointG.array() >= mGN.mSz.array()).any()) {
      std::cout << "Error: basePointG = " << mBasePointG.transpose() << std::endl;
      throw std::runtime_error("Chosen base point is not within the grid");
    }
    l = -mBasePointG;
    u = mGN.mSz - mBasePointG;
    mBasePointIndex = (mBasePointG.array() * mGN.mCumSz.array()).sum();
  }

  void index2grid(int Index) const {
    // Compute grid coordinates from the index
    // note mTmpP is mutable, so const is OK
    for (int i = mTmpP.size()-1; i >= 0; i--) {
      mTmpP[i] = Index/mGN.mCumSz[i];
      Index -= mTmpP[i]*mGN.mCumSz[i];
    }
  }
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // PointServerGrid_h

