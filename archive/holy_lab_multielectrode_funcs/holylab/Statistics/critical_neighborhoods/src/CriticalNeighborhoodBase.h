/**  A base "wrapper" class for critical neighborhoods
 * A CriticalNeighborhood class selects among the core algorithm
 * variants in CriticalNeighborhoodAlgorithms.  These are designed to
 * be very short wrapper classes that let you customize behavior to
 * your particular problem type.
 *
 * This is a base class that handles some of the storage variables.  You may
 * be able to customize the behavior simply by deriving from this class.
 *
 * #include "CriticalNeighborhoodBase.h" <BR>
 *
 */



#ifndef CriticalNeighborhoodBase_h
#define CriticalNeighborhoodBase_h

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
class CriticalNeighborhoodBase
{
public:
// TYPEDEFS AND ENUMS
  typedef typename Accumulator::Scalar Scalar;

// LIFECYCLE
  /** Constructor
   * \param ps a PointServer (e.g., PointServerDistanceDirectC)
   * \param acc an Accumulator (e.g., AccumulatorMoments)
   * \note This class is an "algorithm wrapper" and does not implement any storage.  Both ps and acc are modified by the methods of this class.
   */
  CriticalNeighborhoodBase(PointServer& ps,Accumulator& acc,NeighborhoodStatistics<Scalar> &ns) : mrPS(ps), mrAcc(acc), mrNS(ns), mAccTmp(acc) {;}

  // Default destructor


// OPERATIONS
  /** Inflate the neighborhood until statistical significance is achieved */
  bool inflateToCriterion() {
    return CriticalNeighborhoodAlgorithms::inflateToCriterion(mT2,mSortOrder,mTiesPrevious,mrAcc,mrPS,mrNS);
  }

  /** Inflate over all points */
  void inflateAll() {
    CriticalNeighborhoodAlgorithms::inflateAll(mT2,mSortOrder,mTiesPrevious,mrAcc,mrPS);
  }

  int backtrack1() {
    return CriticalNeighborhoodAlgorithms::backtrackFixedWithTies(mrAcc,mrPS,mSortOrder,mTiesPrevious,mSortOrder.size()-1);
  }

  /** Backtrack to peel off the neighbors that are not part of this neighborhood.  This is an operation that should only be done once you've reached the peak, e.g., after flowToPeak.  On exit, the Accumulator is in a state representing just the "true" neighbors of the peak, and the first n (where n is the return value of this function) elements of sortOrder will list its neighbors.  (In contrast, the PointServer will be in a "corrupted" state, based around a new flow point after visiting the peak.)
   * \return The number of neighbors after backtracking. In general, this will be shorter than sortOrder; only the first n values (where n is the return value) should be considered neighbors.
   */ 
  int backtrackPeak(void)
  {
    // If we exited due to an outlier or exhausting the data set, then
    // don't backtrack
    if (mrPS.isAtEnd())
      return mrAcc.n();
    // Shift to the new centroid, and compute T2 for as many neighbors
    // as were in the neighborhood at the peak. Do all the work in
    // temporary storage, so we don't overwrite the values obtained at
    // the peak.
    mrPS.update(mrAcc);
    CriticalNeighborhoodAlgorithms::inflateAll(mT2Tmp,mSortOrderTmp,mTiesPreviousTmp,mAccTmp,mrPS,mrAcc.n());
    // Find the number of neighbors for which T2 was _most_
    // significant (as estimated by the largest ratio relative to
    // threshold).  Choosing the most-significant brings the centroid
    // "farthest" back in the direction of the peak.
    int nNbrs = 0;
    Scalar W = 0;
    int dof = mrAcc.dof();  // we'll use the "final" dof for all comparisons
    Scalar T2ratio_max = 0;
    mTiesPrevious.push_back(false); // sentinel value; we'll only consider ends of blocks of ties
    for (int n = 0; n < mT2Tmp.size(); n++) {
      W += mrPS.weight(mSortOrderTmp[n]);
      Scalar thisratio = mT2Tmp[n]/mrNS.threshold(dof,std::floor(W+0.5));
      if (thisratio > T2ratio_max && !mTiesPrevious[n+1]) {
	T2ratio_max = thisratio;
	nNbrs = n+1;
      }
    }
    if (T2ratio_max == 0)
      throw std::runtime_error("No suitable backtrack found");
    // Backtrack the peak's accumulator to this number of points
    CriticalNeighborhoodAlgorithms::backtrackFixed(mrAcc,mrPS,mSortOrder,mTiesPrevious,nNbrs);
    mTiesPrevious.pop_back();  // remove the sentinel
    return nNbrs;
  }


// ACCESS
  PointServer& pointServer() { return mrPS; }

  Accumulator& accumulator() { return mrAcc; }

  const std::vector<Scalar>& T2() const { return mT2; }

  const std::vector<int>& sortOrder() const { return mSortOrder; }

  const std::vector<bool>& tiesPrevious() const { return mTiesPrevious; }


// INQUIRY

protected:
private:
  PointServer&                    mrPS;
  Accumulator&                    mrAcc;
  NeighborhoodStatistics<Scalar>& mrNS;
  std::vector<Scalar>             mT2;
  std::vector<int>                mSortOrder;
  std::vector<bool>               mTiesPrevious;
  // Extra storage for backtracking
  Accumulator                     mAccTmp;
  std::vector<Scalar>             mT2Tmp;
  std::vector<int>                mSortOrderTmp;
  std::vector<bool>               mTiesPreviousTmp;

  // Make copy constructor and assignment operator private to prevent
  // copying (at least, until we know we need copying)
  CriticalNeighborhoodBase(const CriticalNeighborhoodBase& from);
  CriticalNeighborhoodBase& operator=(const CriticalNeighborhoodBase& from);
};


#endif  // CriticalNeighborhoodBase_h
