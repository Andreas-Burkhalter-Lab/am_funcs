/**  \brief Computes likelihood that more-distant points are outliers
 *
 * #include "OutlierStatistics.h" <BR>
 *
 * This class evaluates the probability that, given a set of n points
 * drawn from a Gaussian distribution centered at zero, there would
 * likely be yet more points farther from zero.  See T. E. Holy,
 * "Local outlier detection by inﬂation of a bounded domain".
 *
 * This class acts as a cache for evaluating the thresholds needed to
 * achieve statistical significance for any given dof (degrees of
 * freedom) and n (number of data points).  The evaluation of the
 * requisite special functions can take a long time in comparison to
 * other operations, and since the same values may be needed
 * repeatedly, it makes sense to save the results.
 *
 * \tparam Scalar The scalar type used to represent T2. This will
 * presumably be float or double; unsigned types are not supported.
 *
 * \sa NeighborhoodStatistics.h
 */

#ifndef OutlierStatistics_h
#define OutlierStatistics_h

// SYSTEM INCLUDES
//
#include <vector>
#include <map>
#include <boost/math/distributions/gamma.hpp>

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//


template <typename Scalar>
class OutlierStatistics
{
public:
// TYPEDEFS AND ENUMS

// LIFECYCLE

  /** Constructor
   * \param pvalue The p-value that will be used in testing significance
   * \param nmin (optional, default 5) The minimum number of points in a neighborhood before one starts testing for outliers.
   */
  OutlierStatistics(Scalar pvalue, int nmin = 5) : mPValue(pvalue), mNMin(nmin) {;}

  // Prevent copying and assignment by making private

  // Default destructor
  //~OutlierStatistics(void);


// OPERATORS

// OPERATIONS
  /** Find the last non-outlier in a sequence of points
   * \param order, a vector containing the indices of the points in the order that they should be considered
   * \param dist2, a vector of squared-distances to each point (the ith point has distance dist2[i])
   * \param d, the dimensionality
   * \return the number points before reaching the first outlier
   */
  int checkSequence(const std::vector<int>& order,const std::vector<Scalar>& dist2,int d) {
    if (d < 1)
      return order.size();

    Scalar r2,R2cum;
    int n;
    std::vector<Scalar>& v = findTable(d);
    if (order.size() > v.size()) {
      // Fill all missing values with "undefined"
      v.insert(v.end(),order.size()-v.size()+1,undefined);
    }

    R2cum = 0;
    for (n = 1; n <= order.size(); n++) {
      r2 = dist2[order[n-1]];
      if (n >= mNMin) {
	if ( isSignificantCached(v,r2/(R2cum/n),d,n) )
	  break;
      }
      R2cum += r2;
    }
    return n-1;
  }

  // Used for debugging
  // void print() {
  //   typename std::vector< Scalar >::iterator itn;
  //   DofMapType::iterator itdof;
  //   for (itdof = mDofMap.begin(); itdof != mDofMap.end(); itdof++) {
  //     cout << "Dof " << itdof->first << ":" << endl;
  //     for (itn = mNMap[itdof->second].begin(); itn != mNMap[itdof->second].end(); itn++)
  // 	cout << *itn << ' ';
  //     cout << endl;
  //   }
  // }
    

// ACCESS

// INQUIRY
  /** Test whether we expect points beyond a given distance
   * \param dist2 The square Mahalanobis distance to a candidate outlier, where the distance is calculated using the covariance matrix of the points observed so far; alternatively, this is the ratio of the square distance to candidate outlier and the mean square distance of the included points.
   * \param dof The dimensionality
   * \param n The number of points observed so far, deemed to be within the Gaussian
   */
  bool isSignificant(Scalar dist2,int dof,int n) {
    if (n < mNMin || dof < 1)
      return false;
    std::vector<Scalar>& v = findTable(dof);
    return (dist2 > threshold_cached(v,dof,n));
  }

protected:
private:
  enum { undefined = -1};
  typedef typename std::map<int,int> DofMapType;

  Scalar                   mPValue;
  DofMapType               mDofMap;
  std::vector<std::vector<Scalar> >  mNMap;
  int                      mNMin;

  std::vector<Scalar>& findTable(int dof) {
    // Determine whether this value of dof has been seen before
    DofMapType::iterator itdof = mDofMap.find(dof);
    if (itdof == mDofMap.end()) {
      // Haven't seen it, create a new lookup table for it
      mDofMap[dof] = mNMap.size();
      itdof = mDofMap.find(dof);
      mNMap.resize(mNMap.size()+1);
    }
    // Return the corresponding n lookup table
    return mNMap[itdof->second];
  }

  Scalar threshold_cached(std::vector<Scalar>& v,int dof,int n) {
    if (n+1 > v.size()) {
      // This value of n has not been encountered before. Fill all
      // missing values with "undefined"
      v.insert(v.end(),n-v.size()+1,undefined);
    }
    if (v[n] == undefined) {
      // Calculate the threshold for significance
      if (mPValue > 0.05)
	v[n] = boost::math::gamma_p_inv(Scalar(dof)/2,pow(1-mPValue,1.0/(n+1))) * (2.0/dof);
      else {
	// For small pvalue, this version is preferred because it can
	// be more accurate. For pvalue = 0.05 (not terribly small),
	// it seems to be good to 5-6 digits.
	Scalar x = mPValue/(n+1) * (1 + n*mPValue/2/(n+1) * (1 + (2*n+1)*mPValue/3/(n+1) * (1 + (3*n+2)*mPValue/4/(n+1))));  // Taylor expansion of 1 - (1-pvalue)^(1/(n+1))
	v[n] = boost::math::gamma_q_inv(Scalar(dof)/2,x) * (2.0/dof);
      }
    }
    return v[n];
  }

  // Make these private to prevent copying until we know we need this functionality
  OutlierStatistics(const OutlierStatistics& from);
  OutlierStatistics& operator=(const OutlierStatistics& from);  
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // OutlierStatistics_h
