/**  \brief Store the history of point neighborhoods, and test for cycling
 *
 * \details #include "NeighborhoodHistory.h" <BR>
 *
 * Under mean shift, a probe point follows a "trajectory" defined by its neighborhood at each iteration.  This class keeps track of the history of neighbors in the neighborhood, and can indicate whether this trajectory has returned to a previously-visited "peak" of the cycle (defined as the base point with the most neighbors).
 *  
 * @see CriticalNeighborhood
 */

#ifndef NeighborhoodHistory_h
#define NeighborhoodHistory_h

// SYSTEM INCLUDES
//
#include <stdint.h>
#include <vector>

// PROJECT INCLUDES
//

// LOCAL INCLUDES
//

// FORWARD REFERENCES
//


class NeighborhoodHistory
{
public:
// LIFECYCLE
  /** Default constructor. */
  NeighborhoodHistory() {;}
  // Copy constructor made private to prevent copies
  // Default destructor suffices


// OPERATIONS
  /** Clear the history (start fresh) */
  void restart() { mHashValue.clear(); mN.clear(); }

  /** Add a new neighborhood definition
   * \param sortOrder A list of neighbor indices (as a std::vector<int>)
   * \param sz An optional size, if you do not want to use the full length of sortOrder (the default is to use all of sortOrder)
   */
  void add(const std::vector<int>& sortOrder,int sz) {
    mN.push_back(sz);
    uint64_t thisHash = fletcher64((const uint32_t*)(&sortOrder[0]), sz * sizeof(int)/sizeof(uint32_t));
    mHashValue.push_back(thisHash);
    // Formerly the test for cycling was done here, conditional on
    // cycling not having been set true previously.  However, it was
    // discovered that roundoff error can lead to alternations in the
    // point sequence, and so testing just once for cycling was
    // vulnerable to the code getting stuck in an infinite loop. So
    // now it is done whenever isAtMax is called.
  }
  void add(const std::vector<int>& sortOrder) {
    add(sortOrder,sortOrder.size());
  }

// ACCESS
  const std::vector<int>& n_history() const { return mN; }

// INQUIRY
  /** Test whether there is cycling and whether the current
   *  neighborhood contains the largest number of neighbors since
   *  cycling began
   * \return True if current neighborhood is the cycle peak */
  bool isAtMax() const {
    if (mHashValue.empty())
      return false;
    bool atMax = false;
    uint64_t thisHash = mHashValue.back();
    int thisN = mN.back();
    int n_max = thisN;
    for (int i = mHashValue.size()-2; i >= 0; i--) {
      if (mN[i] > n_max)
	n_max = mN[i];
      if (mHashValue[i] == thisHash && n_max == thisN) {
	atMax = true;
	break;
      }
    }
    return atMax;
  }

protected:
private:
  std::vector<uint64_t> mHashValue;
  std::vector<int>      mN;

  // No need for copies, so make private
  NeighborhoodHistory(const NeighborhoodHistory& from);
  NeighborhoodHistory& operator=(const NeighborhoodHistory& from);

  // Hash computation
  uint64_t fletcher64(const uint32_t *data, size_t len)
  {
    uint64_t sum_sequence = 0;
    uint64_t sum = 0;
    const uint32_t *dataEnd = data+len;
    
    for ( ; data < dataEnd; data++) {
      sum = (sum + *data) % 0xFFFFFFFF;
      sum_sequence = (sum_sequence + sum) % 0xFFFFFFFF;
    }
    
    return (sum_sequence << 32) | sum;
  }
};

#endif  // NeighborhoodHistory_h
