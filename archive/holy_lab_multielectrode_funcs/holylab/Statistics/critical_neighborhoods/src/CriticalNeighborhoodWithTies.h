/**  A class for calculating critical neighborhoods in cases of high symmetry, where there are many points "tied" for the same rank.
 *
 * #include "CriticalNeighborhoodWithTies.h" <BR>
 *
 * \sa CriticalNeighborhood
 */



#ifndef CriticalNeighborhoodWithTies_h
#define CriticalNeighborhoodWithTies_h

// SYSTEM INCLUDES
//
#include <vector>

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "CriticalNeighborhoodAlgorithms.h"

// FORWARD REFERENCES
//

template <typename PointServer,typename Accumulator>
class CriticalNeighborhoodWithTies
{
public:
// TYPEDEFS AND ENUMS
  typedef typename Accumulator::Scalar Scalar;

// LIFECYCLE
  /** Constructor
   * \param ps a PointServer (e.g., PointServerDistanceDirectC)
   * \param acc an Accumulator (e.g., AccumulatorMoments)
   * \param maxTies the maximum number of ties that can occur
   * \note This class is an "algorithm wrapper," and the only internal storage is what's needed for tying points.  Both ps and acc are modified by the methods of this class.
   */
  CriticalNeighborhoodWithTies(PointServer& ps,Accumulator& acc,int maxTies) : mrPS(ps), mrAcc(acc), mX(ps.d(),maxTies), mChisq(maxTies) {;}

  // Default destructor


// OPERATIONS
  /** Calculate chisq for neighbhoorhoods of all sizes
   * Points are added to the neighborhood in the order specified by the PointServer.  After adding each point, both neighborhood and pointwise statistics are re-calculated.
   * \param chisqNeighborhood an STL vector storing the neighborhood statistics
   * \param chisqPoint an STL vector storing the pointwise statistics
   * \param sortOrder an STL vector of ints storing the sequence in which points were added (each point addressed by index)
   * \note For efficiency, it is assumed that the PointServer is already at the beginning.  It is up to you to enforce this yourself, perhaps by setting the basePoint or by calling restart().
   * \sa inflateToCriterion
   */
  void inflateAll(std::vector<Scalar>& chisqNeighborhood,
		  std::vector<Scalar>& chisqPoint,
		  std::vector<int>& sortOrder)
  {
    Scalar chisqN,dof;

    mrAcc.restart(mrPS.basePoint());
    chisqNeighborhood.clear();
    chisqPoint.clear();
    sortOrder.clear();
    chisqNeighborhood.reserve(mrPS.N());
    chisqPoint.reserve(mrPS.N());
    sortOrder.reserve(mrPS.N());
    while (!mrPS.isAtEnd()) {
      mX.col(mTieCounter) = mrPS.currentPoint();
      mTieCounter++;
      mrAcc.addPoint(mrPS.currentPoint());
      sortOrder.push_back(mrPS.currentIndex());
      mrPS++;
      if (mrPS.isAtEnd() || !mrPS.isTiesPrevious()) {
	// Process the block of tying points
	mCumNPerGroup.push_back(mTieCounter+mCumNPerGroup.back());
	mrAcc.statisticsCollection(chisqN,mChisq,dof,mX,mTieCounter);
	chisqNeighborhood.insert(chisqNeighborhood.end(),mTieCounter,chisqN);
	chisqPoint.insert(chisqPoint.end(),&mChisq[0],&mChisq[0]+mTieCounter);
      }
    }
  }

  /** A function expanding the neighborhood until the criterion is triggered
   * Points are added to the neighborhood in order.  As the neighborhood grows, the Accumulator keeps track of its statistical properties; when the chisq() value is above the threshold defined by the p-value, neighborhood growth stops.  At the end, the PointServer is in the state left from when the critical neighborhood critierion was triggered, and the Accumulator holds statistics for all added points.  One can "peel off" the "excess" points from the Accumulator using backtrack.
   * \param chisqNeighborhood a vector containing the neighborhood-statistics chisq value for each point in sortOrder (see below)
   * \param chisqPoint a vector containing the pointwise-statistics chisq value for each point in sortOrder
   * \param sortOrder the output ordering of points returned by the PointServer. The size of sortOrder is the number of points needed to trigger the criterion; however, because criterion may require a certain number of "bad" points to push the displacement to statistical significance, you may want to later remove them (see "backtrack").
   * \param pvalue a scalar of the same type as the chisq() return of acc.  This will be turned into a threshold for chisq, which may depend upon the number of degrees of freedom.  A good default choice is pvalue=0.001
   * \param nMin the minimum number of neighboring points to add before checking statistical significance.  In applications like clustering for which data points are i.i.d. draws from a distribution, you probably want to set nMin = ceil(-log(pvalue)), because this is the minimum number of points for which it is unlikely (p < pvalue) that a bootstrap resampling of the full data set will omit all of these neighbors. However, in other applications (e.g., image processing or filtering), the point count may be guaranteed and you may prefer to keep this as 0.  Default value: 0.
   * \return The number of points recommended to be in the critical neighborhood:  if the return value is n, then the first n points in sortOrder should be used.  Note that return values less than nMin are permitted, so if you want to enforce n >= nMin you need to perform this check yourself.  Use "backtrack" to "peel off" excess neighbors from the Accumulator.
   * \note For efficiency, it is assumed that the PointServer is already at the beginning.  It is up to you to enforce this yourself, perhaps by setting the basePoint or by calling restart().
   * \sa backtrack, inflateAll
   */
  int inflateToCriterion(std::vector<Scalar>& chisqNeighborhood,
			 std::vector<Scalar>& chisqPoint,
			 std::vector<int> &sortOrder,
			 Scalar pvalue,
			 int nMin = 0)
  {
    bool satisfied;
    int n,nRemove;
    Scalar chisqN,dof,dofOld,thresh,*pMax;
    
    satisfied = false;
    n = nRemove = 0;
    mrAcc.restart(mrPS.basePoint());
    sortOrder.clear();
    chisqNeighborhood.clear();
    chisqPoint.clear();
    mCumNPerGroup.clear();
    mCumNPerGroup.push_back(0);
    // Surprisingly, calculation of the chi-squared threshhold (thresh)
    // can contribute noticeably to the total amount of compution
    // in this function.  We try to minimize this by caching the
    // threshold for the previous dof.
    dofOld = -1;      // sentinel value to force computation on first check

    mTieCounter = 0;  // k holds the current column of pc to append to
    
    // Add points until the criterion is triggered
    while (!mrPS.isAtEnd()) {
      mX.col(mTieCounter) = mrPS.currentPoint();
      mTieCounter++;
      mrAcc.addPoint(mrPS.currentPoint());
      sortOrder.push_back(mrPS.currentIndex());
      n++;
      mrPS++;
      if (mrPS.isAtEnd() || !mrPS.isTiesPrevious()) {
	// Process the block of tying points
	mCumNPerGroup.push_back(mTieCounter+mCumNPerGroup.back());
	if (n >= nMin) {
	  // Check the criterion
	  mrAcc.statisticsCollection(chisqN,mChisq,dof,mX,mTieCounter);
	  chisqNeighborhood.insert(chisqNeighborhood.end(),mTieCounter,chisqN);
	  chisqPoint.insert(chisqPoint.end(),&mChisq[0],&mChisq[0]+mTieCounter);
	  // If necessary, (re)compute the chi-squared threshold for significance
	  if (dof != dofOld) {
	    thresh = CriticalNeighborhoodAlgorithms::chisqthresh(dof,pvalue);
	    dofOld = dof;
	  }
	  // Test for termination: neighborhood statistics
	  if (chisqN > thresh) {
	    nRemove = ceil(thresh/dof);  // backtrack # required by Cauchy-Schwarz
	    satisfied = true;
	    break;
	  }
	  // Test for termination: pointwise statistics
	  pMax = std::max_element(&mChisq[0],&mChisq[0]+mTieCounter);
	  if (*pMax > thresh) {
	    nRemove = 1;
	    satisfied = true;
	    break;
	  }
	} else {
	  // Insert 0s for chisq
	  chisqNeighborhood.insert(chisqNeighborhood.end(),mTieCounter,0);
	  chisqPoint.insert(chisqPoint.end(),mTieCounter,0);
	}
	// Prepare for next group of tying points
	mTieCounter = 0;
      }
    }
    
    // Compute the return value
    if (satisfied) {
      n -= nRemove;   // adjust # of neighbors for backtracking
      // Ensure that the # of neighbors fits with a tie boundary
      std::vector<int>::iterator it = mCumNPerGroup.end()-1;
      while (n < *it)
	it--;
      return *it;
    } else
      return n;
  }


  /** A version of inflateToCriterion that does not use pointwise statistics */
  int inflateToCriterion(std::vector<Scalar>& chisqNeighborhood,
			 std::vector<int> &sortOrder,
			 Scalar pvalue,
			 int nMin = 0)
  {
    bool satisfied;
    int n,nRemove;
    Scalar chisqN,dof,dofOld,thresh,*pMax;
    
    satisfied = false;
    n = nRemove = 0;
    mrAcc.restart(mrPS.basePoint());
    sortOrder.clear();
    chisqNeighborhood.clear();
    mCumNPerGroup.clear();
    mCumNPerGroup.push_back(0);
    // Surprisingly, calculation of the chi-squared threshhold (thresh)
    // can contribute noticeably to the total amount of compution
    // in this function.  We try to minimize this by caching the
    // threshold for the previous dof.
    dofOld = -1;      // sentinel value to force computation on first check

    mTieCounter = 0;  // k holds the current column of pc to append to
    
    // Add points until the criterion is triggered
    while (!mrPS.isAtEnd()) {
      mTieCounter++;
      mrAcc.addPoint(mrPS.currentPoint());
      sortOrder.push_back(mrPS.currentIndex());
      n++;
      mrPS++;
      if (mrPS.isAtEnd() || !mrPS.isTiesPrevious()) {
	// Process the block of tying points
	mCumNPerGroup.push_back(mTieCounter+mCumNPerGroup.back());
	if (n >= nMin) {
	  // Check the criterion
	  mrAcc.statistics(chisqN,dof);
	  chisqNeighborhood.insert(chisqNeighborhood.end(),mTieCounter,chisqN);
	  // If necessary, (re)compute the chi-squared threshold for significance
	  if (dof != dofOld) {
	    thresh = CriticalNeighborhoodAlgorithms::chisqthresh(dof,pvalue);
	    dofOld = dof;
	  }
	  // Test for termination
	  if (chisqN > thresh) {
	    nRemove = ceil(thresh/dof);  // backtrack # required by Cauchy-Schwarz
	    satisfied = true;
	    break;
	  }
	} else {
	  // Insert 0s for chisq
	  chisqNeighborhood.insert(chisqNeighborhood.end(),mTieCounter,0);
	}
	// Prepare for next group of tying points
	mTieCounter = 0;
      }
    }
    
    // Compute the return value
    if (satisfied) {
      n -= nRemove;   // adjust # of neighbors for backtracking
      // Ensure that the # of neighbors fits with a tie boundary
      std::vector<int>::iterator it = mCumNPerGroup.end()-1;
      while (n < *it)
	it--;
      return *it;
    } else
      return n;
  }


  /** \brief A function that "peels off" the excess neighbors (ones that triggered the critical neighborhood criterion) from the accumulator
   * The "output" is really the accumulator, which is adjusted by this function so that it no longer includes the excess neighbors.  One can, for example, then read out the centroid of the neighborhood from the accumulator.  (Alternatively, one could bypass this function altogether and recalculate the mean of the critical neighborhood simply using sortOrder, but in general this function may require fewer operations.)  The "peeling off" is typically accomplished by subtraction, so do not call this function more than once (or for more than one value of nNbrs) without a fresh call to triggerCriticalNeighborhoodCriterion.
   * \param sortOrder the order of points served by ps (e.g., as calculated by inflateToCriterion)
   * \param nNbrs the number of "real" neighbors (not including excess neighbors) to include in the neighborhood.  Most likely this will come from the return value of inflateToCriterion, although (for example) you may also want to enforce a certain minimum number of neighbors (e.g., to ensure robustness upon bootstrap resampling).
   * \sa inflateToCriterion
   */
  void backtrack(const std::vector<int>& sortOrder,int nNbrs)
  {
    for (int n = sortOrder.size()-1; n >= nNbrs; n--)
      mrAcc.removePoint(mrPS.point(sortOrder[n]));
  }


// ACCESS
  PointServer& pointServer() { return mrPS; }

  Accumulator& accumulator() { return mrAcc; }

// INQUIRY

protected:
private:
  // References to PointServer and Accumulator
  PointServer& mrPS;
  Accumulator& mrAcc;
  // Internal storage for holding ties
  Eigen::Matrix<Scalar,PointServer::PointType::RowsAtCompileTime,Eigen::Dynamic> mX;
  Eigen::Matrix<Scalar,1,Eigen::Dynamic> mChisq;
  std::vector<int> mCumNPerGroup;
  int mTieCounter;

  // Make copy constructor and assignment operator private to prevent
  // copying, at least until we know we need copying
  CriticalNeighborhoodWithTies(const CriticalNeighborhoodWithTies& from);
  CriticalNeighborhoodWithTies& operator=(const CriticalNeighborhoodWithTies& from);
};


#endif  // CriticalNeighborhoodWithTies_h
