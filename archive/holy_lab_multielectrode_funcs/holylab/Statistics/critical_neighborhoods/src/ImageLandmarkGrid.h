/**  A class that manages image landmarks on a regular grid
 *
 * #include "ImageLandmarkGrid.h" <BR>
 *
 * This creates a set of landmarks and assigns each image pixel to its
 * closest landmark. The landmark points are a (small) subset of the
 * pixels on a regular grid.
 *  
 */

#ifndef ImageLandmarkGrid_h
#define ImageLandmarkGrid_h

// SYSTEM INCLUDES
//

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "PointServerImagePV.h"

// FORWARD REFERENCES
//

template <typename CoordP,typename CoordV>
class ImageLandmarkGrid
{
public:
// TYPEDEFS AND ENUMS
  /** The main scalar type */
  typedef typename CoordP::Scalar Scalar;

  /** The scalar data type used to store image data */
  typedef typename CoordV::Scalar ScalarV;

  /** A type representing spatial position in grid coordinates (integers) */
  typedef typename Eigen::Matrix<int,CoordP::RowsAtCompileTime,1> CoordG;

// LIFECYCLE

  /** Default constructor.
   */
  template <typename MetricType>
  ImageLandmarkGrid(const GridNeighbors<CoordP>& imageGN,Scalar landmarkSeparation,const ScalarV *pImage,int n_values,Scalar valueCoefficient,const MetricType &metric = Metric::Euclidean()) {
    init(imageGN,landmarkSeparation,pImage,n_values,valueCoefficient,metric);
  }
    


  // Default destructor is OK
  //~ImageLandmarkGrid(void);


// OPERATORS
// OPERATIONS
  /** The index of the landmark assigned to a given pixel */
  int pixelLandmarkAssignment(int imIndex) const {
    return mPixelLandmarkAssignment(imIndex);
  }

  /** The index of the pixel used for a given landmark */
  int landmarkPixelIndex(int lmIndex) const {
    return mLandmarkPixelIndex[lmIndex];
  }


  /* The coordinates of the pixel used for a given landmark */
  //  MatrixG::ConstColXpr positionG(int lmIndex) const {
  //    return mLandmarkPositionG.col(lmIndex);
  //  }
    
// ACCESS
  /** A constant reference to the vector of landmark assignments */
  const Eigen::Matrix<int,1,Eigen::Dynamic>& pixelLandmarkAssignment() const {
    return mPixelLandmarkAssignment;
  }

  /** A reference to the vector of pixel indices used for all landmarks */
  const std::vector<int>& landmarkPixelIndex() const {
    return mLandmarkPixelIndex;
  }

  /** The separation of landmarks along each axis, in pixels */
  const CoordG& landmarkSeparationG() const { return mLandmarkSeparationG; }

  /** The landmark grid size */
  const CoordG& size() const { return mGridSz; }

// INQUIRY
  /** The number of landmarks */
  int nLandmarks() const { return mLandmarkPixelIndex.size(); }


protected:
private:
  Eigen::Matrix<int,1,Eigen::Dynamic> mPixelLandmarkAssignment;
  std::vector<int>             mLandmarkPixelIndex;
  CoordG                       mLandmarkSeparationG;
  CoordG                       mGridSz;
  

  template <typename MetricType>
  void init(const GridNeighbors<CoordP>& imageGN,Scalar landmarkSeparation,const ScalarV *pImage,int n_values,Scalar valueCoefficient,const MetricType &metric) {
    int n_dimensions = imageGN.mSz.size();
    int i,j;
    
    // Convert the landmarkSeparation into pixel coordinates
    mLandmarkSeparationG.resize(n_dimensions);
    for (i = 0; i < n_dimensions; i++)
      mLandmarkSeparationG(i) = int(round(landmarkSeparation / imageGN.mG2P(i,i)));

    // Check for a quick exit: is every pixel its own landmark?
    if ((mLandmarkSeparationG.array() == 1).all()) {
      mGridSz = imageGN.mSz;
      mPixelLandmarkAssignment.resize(imageGN.N());
      mLandmarkPixelIndex.resize(imageGN.N());
      for (i = 0; i < imageGN.N(); i++)
	mPixelLandmarkAssignment[i] = mLandmarkPixelIndex[i] = i;
      return;
    }
    
    // Define the landmark grid size
    mGridSz.resize(n_dimensions);
    for (i = 0; i < n_dimensions; i++)
      mGridSz(i) = ceil(Scalar(imageGN.mSz(i))/mLandmarkSeparationG(i));
    CoordG gridCumSz(n_dimensions);
    gridCumSz(0) = 1;
    for (i = 1; i < n_dimensions; i++)
      gridCumSz(i) = gridCumSz(i-1)*mGridSz(i-1);
    // We will use a PointServerImagePV to calculate distances. We
    // must ensure that all image points are within the boundaries of
    // a "fake landmark grid," so we also have a version that is 1 larger.
    CoordG szp1(n_dimensions+1);
    szp1.head(n_dimensions).array() = mGridSz.array() + 1;
    szp1(n_dimensions) = 2;  // sentinel value for carry
    
    // Pull out image data at each landmark location
    Eigen::Matrix<ScalarV,Eigen::Dynamic,Eigen::Dynamic> lmValue(n_values,szp1.prod());
    lmValue.setZero();
    CoordG coord(n_dimensions+1);  // +1 sentinel value for carry
    coord.setZero();
    j = 0;
    while (coord(n_dimensions) == 0) {
      if ((coord.head(n_dimensions).array() < mGridSz.array()).all()) {
	// This is a valid grid point, create a landmark for it
	int thisPixelIndex = (coord.head(n_dimensions).array() *
			      mLandmarkSeparationG.array() * 
			      imageGN.mCumSz.array()).sum();
	mLandmarkPixelIndex.push_back(thisPixelIndex);
	// Fetch its image data
	lmValue.col(j) = Eigen::Map<const Eigen::Matrix<ScalarV,Eigen::Dynamic,1> >(pImage+thisPixelIndex*n_values,n_values,1);
      }
      j++;   // invalid pixels are left with 0s for image data
      // Increment with carry
      coord(0)++;
      i = 0;
      while (coord(i) >= szp1(i)) {
	coord(i) = 0;
	i++;
	coord(i)++;
      }
    }

    // Create a PointServer for the grid of landmark pixels
    typename GridNeighbors<CoordP>::TformType lmG2P(n_dimensions,n_dimensions);
    typename GridNeighbors<CoordP>::TformType lmPSkip(n_dimensions,n_dimensions);
    lmPSkip.setZero();
    lmPSkip.diagonal() = mLandmarkSeparationG.template cast<Scalar>();
    lmG2P = imageGN.mG2P * lmPSkip;
    GridNeighbors<CoordP> lmgn(szp1.head(n_dimensions),imageGN.mR,lmG2P);
    PointServerImagePV<CoordP,CoordV> lmps(lmgn,&lmValue(0,0),n_values,metric);
    lmps.valueCoefficient(valueCoefficient);

    // For each image point, find the closest _valid_ landmark grid point
    int nPixels = imageGN.mSz.prod();
    CoordP coordP(n_dimensions);
    coordP.setZero();
    CoordP lmCoord(n_dimensions);
    mPixelLandmarkAssignment.resize(1,nPixels);
    for (i = 0; i < nPixels; i++) {
      lmCoord.noalias() = imageGN.mG2P * coordP;
      lmps.basePoint(ImagePoint<CoordP,CoordV>(lmCoord,Eigen::Map<const CoordV>(pImage+i*n_values,n_values,1)));
      bool isInvalid = (lmps.currentPointG().array() >= mGridSz.array()).any();
      while (!lmps.isAtEnd() && isInvalid) {
	lmps++;
	isInvalid = (lmps.currentPointG().array() >= mGridSz.array()).any();
      }
      if (isInvalid)
	throw std::runtime_error("No valid neighbor was found");
      // Find the index. We have to compute this from scratch, rather
      // than using lmps.currentIndex(), because lmps is created with
      // +1 pixel boundaries in all coordinates.
      mPixelLandmarkAssignment[i] = (lmps.currentPointG().array() * gridCumSz.array()).sum();
      // Increment the coordinates
      j = 0;
      coordP(0)++;
      while (j < n_dimensions-1 && coordP(j) >= imageGN.mSz(j)) {
	coordP(j) = 0;
	j++;
	coordP(j)++;
      }
    }
  }

  // Make copy constructor and assignment operator private until we
  // know we need this functionality
  ImageLandmarkGrid(const ImageLandmarkGrid& from);
  ImageLandmarkGrid& operator=(const ImageLandmarkGrid& from);  
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // ImageLandmarkGrid_h
