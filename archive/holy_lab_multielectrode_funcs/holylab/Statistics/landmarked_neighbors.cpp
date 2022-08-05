#include "landmarked_neighbors.h"
#include <algorithm>
#include <math.h>

template <class Tdata>
Tdata sqrdist(const Tdata *start1,const Tdata *end1,const Tdata *start2)
{
  Tdata ret;
  Tdata dx;

  ret = 0;
  for (; start1 < end1; start1++,start2++) {
    dx = *start1 - *start2;
    ret += dx*dx;
  }
  return ret;
}

template <class Tdata,class Tint>
void landmarked_neighbors<Tdata,Tint>::initialize(const Tdata *yi,const Tdata *xi,landmarkStruct<Tdata,Tint> &lminfoi)
{
  x = xi;
  lminfo = lminfoi;

  //const Tdata *yEnd = y+lminfoi.d;
  Tdata *lmData;
  Tdata closestLandmarkR2,tmpR2;  // variables for square-distance
  Tdata tmpR;  // variable for distance
  int landmarkIndex;
  int i;

  // Copy the "vantage point," to prevent bugs that arise from the
  // user changing the value of y (this is practical because y is
  // small; x and lminfo are large enough that we don't want to make a
  // copy)
  y.clear();
  y.insert(y.end(),yi,yi+lminfo.d);

  landmarks_sdi.clear();
  point_heap.clear();

  // Step #1: determine the distance of each landmark from the current
  // probe point.  While doing this, keep track of the closest one, so
  // we can return its index
  closestLandmarkR2 = -1;
  //mexPrintf("lminfo.landmarks[0] %g, n_landmarks %d, d %d\n",*(lminfo.landmarks),lminfo.n_landmarks,lminfo.d);
  //mexPrintf("y %x, yEnd %x, diff %d\n",y,y.end(),y.end()-y);

  //mexPrintf("Landmark groups & R2s:\n");
  for (landmarkIndex = 0, lmData = lminfo.landmarks; landmarkIndex < lminfo.n_landmarks; landmarkIndex++, lmData+=lminfo.d) {
    tmpR2 = sqrdist(&(*y.begin()),&(*y.end()),lmData);
    //mexPrintf("%d %g\n",landmarkIndex,tmpR2);
    landmarks_sdi.push_back(SqrdistIndex<Tdata>(tmpR2,landmarkIndex));
    if (closestLandmarkR2 < 0 || tmpR2 < closestLandmarkR2) {
      closestLandmarkR2 = tmpR2;
      closestLandmarkIndex = landmarkIndex;
    }
  }

  if (use_landmarks) {
    // Step #2: For each landmark group, determine (using the triangle
    // inequality) the lower bound on distance between points in the
    // grouping and the probe point.  This is the distance we'll store
    // and sort by.
    for (landmarkIndex = 0; landmarkIndex < lminfo.n_landmarks; landmarkIndex++) {
      tmpR = lminfo.landmarkR[landmarkIndex];  // the radius of the group
      if (landmarks_sdi[landmarkIndex].R2 < tmpR*tmpR) {
	// The current probe point could be inside this landmark
	// group, so the lower bound on distance is 0
	landmarks_sdi[landmarkIndex].R2 = 0;
      } else {
	// The landmark group is farther away, use the triangle
	// inequality to determine the lower bound on the (square)
	// distance
	tmpR = sqrt(landmarks_sdi[landmarkIndex].R2) - tmpR;
	landmarks_sdi[landmarkIndex].R2 = tmpR*tmpR;
      }
    }

    // Step #3: Sort groups by this lower-bound distance (in
    // increasing order).
    //mexPrintf("Unsorted landmark order:\n");
    //for (landmarkIndex = 0; landmarkIndex < lminfo.n_landmarks; landmarkIndex++)
    //  mexPrintf("%d (%g) ",landmarks_sdi[landmarkIndex].index,landmarks_sdi[landmarkIndex].R2);
    //mexPrintf("\n");
    std::sort(landmarks_sdi.begin(),landmarks_sdi.end());
    //mexPrintf("Sorted landmark order:\n");
    //for (landmarkIndex = 0; landmarkIndex < lminfo.n_landmarks; landmarkIndex++)
    //  mexPrintf("%d (%g) ",landmarks_sdi[landmarkIndex].index,landmarks_sdi[landmarkIndex].R2);
    //mexPrintf("\n");

    // Step #4: Clear the heap and initialize the progress variables
    point_heap.clear();
    point_heap_end = point_heap.begin();
    lmIterator = landmarks_sdi.begin();
    currentR2 = lmIterator->R2;

    // Step #5: Get the point_heap set up for the first requests
    advance_heap();
  }


  else {
    // We're not going to be using the landmarks to triage the data
    // points. So instead, compute the distance to each data point and
    // then sort.
    const Tdata *dataEnd = x + lminfo.d * N;
    const Tdata *y_begin = &(*y.begin());
    const Tdata *y_end = &(*y.end());
    const Tdata *dataIterator;
    int dataIndex;
    for (dataIndex = 0, dataIterator = x; dataIterator < dataEnd; dataIterator += lminfo.d, dataIndex++) {
      tmpR2 = sqrdist(y_begin,y_end,dataIterator);
      point_heap.push_back(SqrdistIndex<Tdata>(tmpR2,dataIndex));
    }
    
    std::sort(point_heap.begin(),point_heap.end());

    cur_position = point_heap.begin();
  }
}

template <class Tdata,class Tint>
void landmarked_neighbors<Tdata,Tint>::advance_heap()
{
  int thisGroup,pointIndex,thisPointIndex;
  const Tdata *thisPoint;
  Tdata tmpR2;

  if (!use_landmarks) {
    return;
  }

  // Add points to heap (if necessary).  Here's the scheme &
  // considerations:
  // (1) add all groups that have a lower-bound distance
  //     within the currentR2
  // (2) if the closest point yet has a square distance larger than
  //     currentR2, it might not be the closest point: there might be
  //     another landmark group containing a closer point.  So we increase
  //     currentR2 and add all the groups up to that distance.
  //
  // This leads to a double-while construction.  We could do it with a
  // single-while, but it could result in unnecessary points being
  // added to the heap, if a later group already scheduled for
  // addition has the closest point.
  //mexPrintf("1: currentR2 %g, lmIterator->R2 %g\n",currentR2,lmIterator->R2);
  while (lmIterator < landmarks_sdi.end() && lmIterator->R2 <= currentR2) {
    //mexPrintf("2: currentR2 %g, lmIterator->R2 %g\n",currentR2,lmIterator->R2);
    while (lmIterator < landmarks_sdi.end() && lmIterator->R2 <= currentR2) {
      // Add all the points in the current group
      thisGroup = lmIterator->index;
      //mexPrintf("Adding group %d\n",thisGroup);
      for (pointIndex = 0; pointIndex < lminfo.n_landmarkList[thisGroup]; pointIndex++) {
	thisPointIndex = (int) lminfo.landmarkList[thisGroup][pointIndex] - lminfo.index_offset;
	thisPoint = x + lminfo.d*thisPointIndex;
	tmpR2 = sqrdist(&(*y.begin()),&(*y.end()),thisPoint);
	*point_heap_end = SqrdistIndex<Tdata>(tmpR2,thisPointIndex);
	// make it so the top of the heap is the closest point
	push_heap(point_heap.begin(),++point_heap_end,is_farther<Tdata>());
      }
      //mexPrintf("Heap status:\n");
      //typename vector< SqrdistIndex<Tdata> >::iterator phI;
      //for (phI = point_heap.begin(); phI < point_heap_end; phI++)
      //  mexPrintf("x%d (%g) ",phI->index,phI->R2);
      //mexPrintf("\n");
      lmIterator++;
    } // end of inside while loop
    // We want to make sure we keep including groups until we get some
    // points on the heap that are closer than any points not on the
    // heap.  The way to do that is to adjust currentR2 to reflect the
    // distance to the closest point.
    if (!is_empty())
      currentR2 = point_heap.begin()->R2;
    else if (lmIterator < landmarks_sdi.end()) {
      // If the heap is still empty, but we haven't exhausted all
      // landmarks, then we must have had only empty landmark
      // groups. In that case, increment to the next landmarkR2.
      currentR2 = lmIterator->R2;
    }
  }
}

template <class Tdata,class Tint>
void landmarked_neighbors<Tdata,Tint>::next_point()
{
  if (!use_landmarks) {
    cur_position++;
    return;
  }

  pop_heap(point_heap.begin(),point_heap_end,is_farther<Tdata>());
  point_heap_end--;
 
  if (!is_empty())
    currentR2 = point_heap.begin()->R2;
  else if (lmIterator < landmarks_sdi.end())
    currentR2 = lmIterator->R2;
  advance_heap();
}

