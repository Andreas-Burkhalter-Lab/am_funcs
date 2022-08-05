/**  \brief Computes statistical significance of displacement of the mean
 *
 * #include "NeighborhoodStatistics.h" <BR>
 *
 * This class evaluates the signficance of the displacement of the
 * mean.  In general this is Hotelling's T2 statistic, but with
 * constrained covariance models (e.g., Isotropic or Diagonal) we need
 * to use the appropriate statistic.  This implements the results in
 * T. E. Holy, "Statistics of the mean using constrained sample
 * covariance matrices".
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
 * \note For diagonal covariances, grouping is not yet supported.
 *  
 * \sa OutlierStatistics.h
 */

#ifndef NeighborhoodStatistics_h
#define NeighborhoodStatistics_h

// SYSTEM INCLUDES
//
#include <vector>
#include <map>
#include <boost/math/distributions/fisher_f.hpp>

// PROJECT INCLUDES
//
#include "AccumulatorMoments.h"

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//


template <typename Scalar>
class NeighborhoodStatistics
{
public:
// TYPEDEFS AND ENUMS

// LIFECYCLE

  /** Constructor
   * \param pvalue The p-value that will be used in testing significance.  Should be a scalar of the same type as the T2 return of acc.  This will be turned into a threshold for T2, which may depend upon the number of degrees of freedom.  A good default choice is pvalue=0.001.
   * \param nMin the minimum number of neighboring points to add before checking statistical significance.  In applications like clustering for which data points are i.i.d. draws from a distribution, you probably want to set nMin = ceil(-log(pvalue)), because this is the minimum number of points for which it is unlikely (p < pvalue) that a bootstrap resampling of the full data set will omit all of these neighbors. However, in other applications (e.g., image processing or filtering), the point count may be guaranteed and you may prefer not to constrain this value.  Default value: 2 (and any smaller value will be set to 2, by necessity for the division by n-1).
   * \param ctype A flag indicating the covariance model to use in testing the significance of the displacement
   */
  NeighborhoodStatistics(Scalar pvalue,int nmin,Moments::CovarianceModel ctype = Moments::Isotropic) : mPValue(pvalue), mNMin(nmin), mCovarianceModel(ctype) {
    set_max();
    if (mNMin < 2)
      mNMin = 2;
  }

  NeighborhoodStatistics(Scalar pvalue,Moments::CovarianceModel ctype = Moments::Isotropic) : mPValue(pvalue), mCovarianceModel(ctype) {
    set_max();
    mNMin = 2;
  }

  // Prevent copying and assignment by making private

  // Default destructor
  //~NeighborhoodStatistics(void);


// OPERATORS

// OPERATIONS
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
  Scalar threshold(int dof,int n) {
    if (n < mNMin || dof < 1)
      return mMax;
    // Get the corresponding n lookup table
    std::vector<Scalar>& v = findTable(dof);
    // Compare the threshold
    return threshold_cached(v,dof,n);
  }

// INQUIRY
  /** Test whether the statistic is significant
   * \param T2 Hotelling's T2 statistic for displacement of the mean
   * \param dof The dimensionality
   * \param n The number of points in the neighborhood
   */
  bool isSignificant(Scalar T2,int dof,int n) {
    if (n < mNMin || dof < 1)
      return false;
    // Get the corresponding n lookup table
    std::vector<Scalar>& v = findTable(dof);
    // Compare the threshold
    return (T2 > threshold_cached(v,dof,n));
  }


protected:
private:
  enum { undefined = -1};
  typedef typename std::map<int,int> DofMapType;

  Scalar                   mPValue;
  Moments::CovarianceModel mCovarianceModel;
  DofMapType               mDofMap;
  std::vector<std::vector<Scalar> >  mNMap;
  Scalar                   mMax;
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
      Scalar s,n2;
      if (mCovarianceModel == Moments::Isotropic) {
	s = Scalar(n-1)/(dof*n-1);
	n2 = dof*(n-1);
      } else if (mCovarianceModel == Moments::Diagonal) {
	// Diagonal covariance. Here we use the 2-moment approximation
	if (n < 5) {
	  v[n] = mMax;
	  return false;
	} else {
	  n2 = Scalar((dof+2)*n-5*dof+2)/3;
	  s = Scalar((n-3)*n2)/(dof*(n-1)*(n2-2));
	}
      } else {
	// Full covariance model
	if (n <= dof) {
	  v[n] = mMax;
	  return false;
	} else {
	  s = Scalar(n-dof)/(dof*(n-1));
	  n2 = n-dof;
	}
      }
      v[n] = boost::math::quantile(boost::math::fisher_f(dof,n2),1-mPValue)/s;
    }
    return v[n];
  }

  void set_max() {
    if (std::numeric_limits<Scalar>::has_infinity)
      mMax = std::numeric_limits<Scalar>::infinity();
    else
      mMax = std::numeric_limits<Scalar>::max();
  }

  // Make these private to prevent copying until we know we need this functionality
  NeighborhoodStatistics(const NeighborhoodStatistics& from);
  NeighborhoodStatistics& operator=(const NeighborhoodStatistics& from);  
};

// INLINE METHODS
//

// EXTERNAL REFERENCES
//

#endif  // NeighborhoodStatistics_h
