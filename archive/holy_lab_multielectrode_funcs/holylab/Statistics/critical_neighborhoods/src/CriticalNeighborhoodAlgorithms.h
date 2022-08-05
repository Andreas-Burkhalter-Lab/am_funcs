#ifndef CriticalNeighborhoodAlgorithms_h
#define CriticalNeighborhoodAlgorithms_h

// SYSTEM INCLUDES
//
#include <vector>
#include <stdexcept>
#include <boost/math/distributions/gamma.hpp>

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//
#include "NeighborhoodHistory.h"
#include "NeighborhoodStatistics.h"

// FORWARD REFERENCES
//

/// \namespace CriticalNeighborhoodAlgorithms 
/// \brief Define and manipulate critical neighborhoods
/// \details This header file defines generic functions needed to calculate the "critical neighborhood" for a point.  A neighborhood grows by accumulating points, for example in order of increasing distance from the base point.  A "critical" neighborhood is one that has only just triggered some statistical criterion for showing a trend. \n The ordering of points is determined by a PointServer, and the statistical properties of the region are measured by the Accumulator.  The Accumulator returns T2() and the degrees of freedom dof(), which can be used to judge statistical significance.

namespace CriticalNeighborhoodAlgorithms
{
  /** Calculate T2 for neighbhoorhoods of all sizes
   * Points are added to the neighborhood in the order specified by the PointServer.  After adding each point, both neighborhood and pointwise statistics are re-calculated.
   * \param T2 an STL vector storing the neighborhood statistics
   * \param sortOrder an STL vector of ints storing the sequence in which points were added (each point addressed by index)
   * \param tiesPrevious a vector, tiesPrevious[i]=true indicates that point sortOrder[i] should be considered as having the same rank as point sortOrder[i-1].
   * \param acc an Accumulator (e.g., AccumulatorMoments)
   * \param ps a PointServer (e.g., PointServerDistanceDirectC)
   * \param nMax the maximum number of neighbors to include. -1 (the default) indicates that all neighbors should be used.
   * \note If your points are weighted, then any statistical test you perform on T2[i] should be using n = the sum of the weights for points sortOrder[0...i], rather than n = i.  Likewise, nMax constrains the total sum of the weights, not the number of points.
   * \note For efficiency, it is assumed that the PointServer is already at the beginning.  It is up to you to enforce this yourself, perhaps by setting the basePoint or by calling restart().
   * \sa inflateToCriterion
   */
  template <typename AccumulatorType,typename PointServerType>
  void inflateAll(std::vector<typename AccumulatorType::Scalar>& T2,
		  std::vector<int>& sortOrder,
		  std::vector<bool> &tiesPrevious,
		  AccumulatorType &acc,
		  PointServerType &ps,
		  typename AccumulatorType::Scalar nMax = -1)
  {
    typename AccumulatorType::Scalar tmpT2,dof;
    int i,N;
    int firstInGroup;  // first element in a sequence of tying points

    // Initialize
    acc.restart(ps.basePoint());
    T2.clear();
    sortOrder.clear();
    tiesPrevious.clear();
    N = ps.N();
    //    Commented out because we want the sum of weights to be nMax,
    //    not the number of entries to be nMax
    //    if (nMax > 0)
    //      N = (nMax < N) ? nMax : N;
    T2.reserve(N);
    sortOrder.reserve(N);
    tiesPrevious.reserve(N);
    firstInGroup = 0;

    // Add points to accumulator in blocks of tying points
    while (!ps.isAtEnd()) {
      if (nMax > 0 && acc.n() >= nMax)
	break;
      // Store the sequence
      sortOrder.push_back(ps.currentIndex());
      tiesPrevious.push_back(ps.isTiesPrevious());
      // Add the point to the accumulator
      ps.addCurrentPoint(acc);
      // Increment the PointServer
      ps++;
      if (ps.isAtEnd() || !ps.isTiesPrevious() || acc.n() >= nMax) {
	// Compute neighborhood statistics
	acc.statisticsNeighborhood(tmpT2,dof);
	T2.insert(T2.end(),sortOrder.size()-firstInGroup,tmpT2);
	firstInGroup = sortOrder.size();  // reset for next tying group
      }
    }
  }

  /** A function expanding the neighborhood until the criterion is triggered
   * Points are added to the neighborhood in order.  As the neighborhood grows, the Accumulator keeps track of its statistical properties; when the T2() value is above the threshold defined by the p-value, neighborhood growth stops.  At the end, the PointServer is in the state left from when the critical neighborhood critierion was triggered, and the Accumulator holds statistics for all added points.  One can "peel off" the "excess" points from the Accumulator using backtrack.
   * \param T2 a vector containing the neighborhood-statistics T2 value for each point in sortOrder (see below)
   * \param sortOrder the output ordering of points returned by the PointServer. The size of sortOrder is the number of points needed to trigger the criterion; however, because criterion may require a certain number of "bad" points to push the displacement to statistical significance, you may want to later remove them (see "backtrack").
   * \param tiesPrevious a vector, tiesPrevious[i]=true indicates that point sortOrder[i] should be considered as having the same rank as point sortOrder[i-1].
   * \param acc the Accumulator (e.g., AccumulatorMoments)
   * \param ps the PointSever (e.g., PointServerDistanceDirectC)
   * \return The flag indicating the reason that no more points were added to the neighborhood: 0 = all data points added; 1 = the PointServer indicated that an outlier had been detected; 2 = the neighborhood criterion was met.  The accumulator is left in a state including all tested points (to modify the accumulator, see "backtrack" functions).
   * \note If your points are weighted, then the significance of T2[i] is evaluated using n = the sum of the weights for points sortOrder[0...i], rather than n = i.
   * \note For efficiency, it is assumed that the PointServer is already at the beginning.  It is up to you to enforce this yourself, perhaps by setting the basePoint or by calling restart().
   * \sa backtrack, inflateAll
   */
  template <typename AccumulatorType,typename PointServerType>
  bool inflateToCriterion(std::vector<typename AccumulatorType::Scalar>& T2,
			  std::vector<int> &sortOrder,
			  std::vector<bool> &tiesPrevious,
			  AccumulatorType &acc,
			  PointServerType &ps,
			  NeighborhoodStatistics<typename AccumulatorType::Scalar> &nbrhStats)
  {
    bool satisfied;
    int n;
    int firstInGroup;  // first element in a sequence of tying points
    typename AccumulatorType::Scalar tmpT2,dof;
    
    satisfied = false;
    acc.restart(ps.basePoint());
    T2.clear();
    sortOrder.clear();
    tiesPrevious.clear();
    firstInGroup = 0;
    
    // Add points until the criterion is triggered
    while (!ps.isAtEnd()) {
      ps.addCurrentPoint(acc);
      sortOrder.push_back(ps.currentIndex());
      tiesPrevious.push_back(ps.isTiesPrevious());
      ps++;
      if (ps.isAtEnd() || !ps.isTiesPrevious()) {
	// Compute neighborhood statistics
	acc.statisticsNeighborhood(tmpT2,dof);
	if (!std::isfinite(tmpT2))
	  throw std::runtime_error("tmpT2 is not finite");
	T2.insert(T2.end(),sortOrder.size()-firstInGroup,tmpT2);
	if (nbrhStats.isSignificant(tmpT2,dof,std::floor(acc.n()+0.5))) {
	  satisfied = true;
	  break;
	}
	// Prepare for the next tying group
	firstInGroup = sortOrder.size();
      }
    }
    
    return satisfied;
  }


  /** \brief A function that "peels off" the excess neighbors (ones that triggered the critical neighborhood criterion) from the accumulator
   * The "output" is the accumulator, which is adjusted by this function so that it no longer includes the excess neighbors.  One can, for example, then read out the centroid of the neighborhood from the accumulator.  (Alternatively, one could bypass this function altogether and recalculate the mean of the critical neighborhood simply using sortOrder, but in general this function may require fewer operations.)  The "peeling off" is typically accomplished by subtraction, so if you call it twice (without altering sortOrder), you will double-count some removed points.
   * \param acc the Accumulator (e.g., AccumulatorMoments)
   * \param ps the PointSever (e.g., PointServerDistanceDirectC)
   * \param sortOrder the order of points served by ps (e.g., as calculated by inflateToCriterion)
   * \param tiesPrevious the same-named output of inflateToCriterion. This version ignores information about ties (it's included as an argument for reasons of API compatbility), but see backtrackFixedWithTies for a backtracking algorithm that makes used of this.
   * \param nNbrs the number of "real" neighbors (not including excess neighbors) to include in the neighborhood.  Most likely this will come from the return value of inflateToCriterion, although (for example) you may also want to enforce a certain minimum number of neighbors (e.g., to ensure robustness upon bootstrap resampling).
   * \sa inflateToCriterion,backtrackFixedWithTies
   */
  template <typename AccumulatorType,typename PointServerType>
  int backtrackFixed(AccumulatorType& acc,const PointServerType& ps,const std::vector<int>& sortOrder,const std::vector<bool>& tiesPrevious,int nNbrs)
  {
    int n;
    for (n = sortOrder.size()-1; n >= nNbrs; n--)
      acc.removePoint(ps.point(sortOrder[n]));
    return n+1;  // +1 because for loop terminates when n = nNbrs-1
  }

  /** \brief A version of backtrackFixed that also removes any additional points in the same "tying" group
   * \sa backtrackFixed,inflateToCriterion
   */
  template <typename AccumulatorType,typename PointServerType>
  int backtrackFixedWithTies(AccumulatorType& acc,const PointServerType& ps,const std::vector<int>& sortOrder,const std::vector<bool>& tiesPrevious,int nNbrs)
  {
    int n;
    for (n = sortOrder.size()-1; n >= nNbrs; n--) {
      acc.removePoint(ps.point(sortOrder[n]));
    }
    while (n > 0 && tiesPrevious[n]) {
      acc.removePoint(ps.point(sortOrder[n]));
      n--;
    }
    return n+1; // +1 because for loop terminates when n is not removed
  }

  /** \brief Compute "correlations" between individual data points and the mean displacement from the centroid.
   * This will someday be used for a backtrack algorithm that finds the most "problematic" points and eliminates them.
   * \param corr a vector holding the "correlation" (really the dot product between the point and the mean displacement). Note that the dot product is with the point, and not the point's displacement from the base point. (It is possible that this might change in the future, if roundoff is a problem.)
   * \param temp storage for a single displacement
   * \param acc the Accumulator (e.g., AccumulatorMoments)
   * \param ps the PointSever (e.g., PointServerDistanceDirectC)
   * \param sortOrder the order of points served by ps (e.g., as calculated by inflateToCriterion)
   */
  template <typename AccumulatorType,typename PointServerType>
  void computeCorrelations(std::vector<typename AccumulatorType::Scalar>& corr,typename PointServerType::PointType& temp,AccumulatorType& acc,const PointServerType& ps,const std::vector<int>& sortOrder)
  {
    temp = acc.mean() - acc.basePoint();
    acc.solveInPlace(temp);
    corr.resize(sortOrder.size());
    for (int i = 0; i < sortOrder.size(); i++)
      corr[i] = temp.dot(ps.point(sortOrder[i]));
  }

  // FIXME create backtrackCorrelated, removing the given # of points that are most highly correlated with the mean. For this, you want a "solve" method in the accumulator so that inversion only has to be done once (i.e., don't use dotProduct directly, just use dot() on a "solved" mean vector)



  /** \brief Flow a point to its peak
   * Given a starting point, this algorithm iteratively performs mean shift until the peak is reached.
   * \param cn the CriticalNeighborhood structure, which selects the variants of algorithms (inflateToCriterion and backtrack) you wish to use for each step of flow. See CriticalNeighborhoodBase for a base class.
   * \param history a class (and workspace) for testing cycling
   * \return The number of points in the neighborhood
   */
  template <typename CNType>
  int flowToPeak(CNType& cn,NeighborhoodHistory& history)
  {
    int nNbrs;
    bool isAtMax = false;

    history.restart();
    while (true) {
      // Compute the critical neighborhood
      cn.inflateToCriterion();
      nNbrs = cn.sortOrder().size();
      //nNbrs = cn.backtrack1();
      // Determine whether we are at the peak
      history.add(cn.sortOrder(),nNbrs);
      isAtMax = history.isAtMax();
      if (isAtMax)
	break;
      // Shift the base point (and any other updating)
      cn.pointServer().update(cn.accumulator());
    }
    return nNbrs;
  }
    
  struct FlowResults {
    int nNbrs;
    bool isAtMax;
  };

  /** \brief Flow a point to its peak, or until the nearest-neighbor changes
   * Given a starting point, this algorithm iteratively performs mean shift until the peak is reached or the nearest-neighbor changes.
   * \param cn the CriticalNeighborhood structure, which selects the variants of algorithms (inflateToCriterion and backtrack) you wish to use for each step of flow. See CriticalNeighborhoodBase for a base class.
   * \param history a class (and workspace) for testing cycling
   * \param pvalue the statistical significance required to reach the critical neighborhood criterion
   * \param nMin the minimum number of points in the critical neighborhood (default 0)
  * \return a FlowResults structure, holding the number (nNbrs) of neighbors at termination and a flag (isAtMax) reporting whether termination was caused by reaching the peak
*/
  template <typename CNType>
  FlowResults flowToNeighbor(CNType& cn,NeighborhoodHistory& history,typename CNType::Scalar pvalue,int nMin = 1)
  {
    int nNbrs,firstIndex;
    bool isAtMax = false;
    bool isFirst = true;

    history.restart();
    while (true) {
      // Compute the critical neighborhood
      nNbrs = cn.inflateToCriterion(pvalue,nMin);
      if (nNbrs < nMin)
	nNbrs = nMin;
      cn.backtrack(nNbrs);
      // If this is the first iteration, store the nearest-neighbor index
      if (isFirst) {
	firstIndex = cn.sortOrder()[0];
	isFirst = false;
      }
      // Determine whether the nearest-neighbor has changed
      if (cn.sortOrder()[0] != firstIndex)
	break;
      // Determine whether we are at the peak
      history.add(cn.sortOrder(),nNbrs);
      isAtMax = history.isAtMax();
      if (isAtMax)
	break;
      // Shift the base point (and perform any other updating)
      cn.pointServer().update(cn.accumulator());
    }

    FlowResults results;
    results.nNbrs = nNbrs;
    results.isAtMax = isAtMax;
    return results;
  }
  

  /// \brief A function to assign data points to peaks
  /// Under mean shift, probe points move to a peak. Rather than flowing each point until it stops moving, it can be much more efficient to terminate flow once it can be associated with another probe point, so that the two have a common "destiny."  This results in a "map" of neighbor associations; this function "flows the map" to determine the final "target" of each probe point.
  /// \param map the initial associations as a vector of indices. map[i] = j implies that the ith point is to be associated with the jth one; j is "uphill" of i.  This parameter is modified; on output, it contains the final map of associations, i pointing to the probe point closest to the "peak."
  /// \param n a vector containing the number of points in the neighborhood for each point. This is used to break cycles and/or ties, to help determine the true "peaks."
  void flowMap(std::vector<int>& mapsTo,const std::vector<int>& n)
  {
    // Copy the input, for later cycle-breaking
    std::vector<int> mapsTo0 = mapsTo;

    // Storage for the unique "targets"
    std::vector<int> tmpMapsTo = mapsTo;
    std::vector<int> tmpMapsToOld(mapsTo.size());
    std::vector<int>::iterator umapsToEnd;
    int uSize,uOldSize;
    std::pair<std::vector<int>::iterator,std::vector<int>::iterator> mismatchPair;

    // Initialize the unique targets
    std::sort(tmpMapsTo.begin(),tmpMapsTo.end());
    umapsToEnd = std::unique(tmpMapsTo.begin(),tmpMapsTo.end());
    uSize = umapsToEnd-tmpMapsTo.begin();

    // Iteratively flow until the map stops changing
    while (true) {
      std::copy(tmpMapsTo.begin(),umapsToEnd,tmpMapsToOld.begin());
      uOldSize = uSize;
      std::copy(mapsTo.begin(),mapsTo.end(),tmpMapsTo.begin());
      std::sort(tmpMapsTo.begin(),tmpMapsTo.end());
      umapsToEnd = std::unique(tmpMapsTo.begin(),tmpMapsTo.end());
      uSize = umapsToEnd-tmpMapsTo.begin();
      if (uSize == uOldSize) {
	// The sizes fit, let's check each value
	mismatchPair = std::mismatch(tmpMapsTo.begin(),umapsToEnd,tmpMapsToOld.begin());
	if (mismatchPair.first == umapsToEnd)
	  break;
      }
    }

    // Because of the possibility of cycles, now we have to repair the map
    // Example: if we have an initial map [1 2 1], i.e.,
    //   0->1, 1->2, 2->1
    // then when iterated this yields [2 1 2], which is stable under
    // further iterations.  However, the final map is disconnected,
    // whereas the original is connected.  We can reconnect by
    // reference to the original map.
    int j,m1,m2,n1,n2;
    int changed = 1;
    while (changed > 0) {
      changed = 0;
      for (int i = 0; i < mapsTo.size(); i++) {
	m1 = mapsTo[i];
	j = mapsTo0[i];
	m2 = mapsTo[j];
	if (m1 != m2) {
	  n1 = n[m1];
	  n2 = n[m2];
	  // Choose the target with the largest neighborhood
	  if (n2 > n1)
	    mapsTo[i] = m2;
	  else if (n2 == n1) {
	    mapsTo[i] = (m1 < m2) ? m1 : m2; // break ties, choose lowest index
	    mapsTo[j] = mapsTo[i];
	  }
	  else
	    mapsTo[j] = m1;
	  changed++;
	}
      }
    }
  }	  
      
  
} // end namespace CriticalNeighborhoodAlgorithms

#endif  // CriticalNeighborhoodAlgorithms_h
