#include <vector>
#include <iterator>

using namespace std;

// A structure to hold landmarks. This is implemented in terms of basic
// pointers so that it's easy to bind to Matlab or other sources
template <class Tdata,class Tint>
struct landmarkStruct {
  int n_landmarks;
  int d;   // The dimensionality
  Tdata *landmarks;   // landmark positions (d-by-n_landmarks)
  Tint *landmarkAssignment; // landmark index corresponding to each point
  Tint **landmarkList;  // List of points associated with each landmark
  int *n_landmarkList;  // Number of points associated with each landmark
  double *landmarkR;   // Radius of each landmark group around the landmark
  int index_offset;   // (For Matlab this should be 1 because it's unit-offset)

  landmarkStruct() {
    n_landmarks = 0;
    d = 0;
    landmarks = NULL;
    landmarkAssignment = NULL;
    landmarkR = NULL;
    landmarkList = NULL;
    n_landmarkList = NULL;
    index_offset = 0;
  }
};



// Class for holding a square-distance and a corresponding index; this
// will be heavily used in heapsort routines
template <class T>
struct SqrdistIndex {
  T R2;   // distance^2
  int index;

  SqrdistIndex() {;}
  SqrdistIndex(T r2,int i) {R2 = r2; index = i;}
  bool operator<(const SqrdistIndex &sdi) const {return (R2 < sdi.R2);}
  bool operator>(const SqrdistIndex &sdi) const {return (R2 > sdi.R2);}
};

template <class T>
struct is_farther : public binary_function<SqrdistIndex<T>,SqrdistIndex<T>,bool> {
  bool operator()(SqrdistIndex<T> &sdi1,SqrdistIndex<T> &sdi2) { return sdi1.R2 > sdi2.R2;}
};


// Class for returning neighbors, in order of distance.  If
// use_landmarks is true, it makes use of landmarked groups of points,
// so that one doesn't have to examine each point.
const bool use_landmarks = 1; // use a const so checks are optimized away
template <class Tdata,class Tint>
class landmarked_neighbors {
 private:
  vector< SqrdistIndex<Tdata> > landmarks_sdi;
  typename vector< SqrdistIndex<Tdata> >::iterator lmIterator;
  int closestLandmarkIndex;
  vector< SqrdistIndex<Tdata> > point_heap;
  typename vector< SqrdistIndex<Tdata> >::iterator point_heap_end;
  Tdata currentR2;  // the radius-squared up to which we're delivering points

  vector<Tdata> y;
  const Tdata *x;
  landmarkStruct<Tdata,Tint> lminfo;
  int N;

  void advance_heap();

  // These are variables that are used when we turn off the triaging
  // by landmark, which we only want to do when we want to find out
  // how much time that's saving us.
  typename vector< SqrdistIndex<Tdata> >:: iterator cur_position;

 public:
  landmarked_neighbors() {;}
  
  // "allocate" has to be called explicitly, before anything else is done
  void allocate(int n_landmarks,int n_points) {
    landmarks_sdi.reserve(n_landmarks);
    point_heap.reserve(n_points);
    N = n_points;
  }
  // "initialize" sets up the data structure that allows one to ask
  // for information about points and landmarks
  void initialize(const Tdata *yi,const Tdata *xi,landmarkStruct<Tdata,Tint> &lminfoi);

  // Information about landmarks
  int closest_landmark() {return closestLandmarkIndex;}
  int landmark_of_current_point() {return lminfo.landmarkAssignment[current_point()->index] - lminfo.index_offset;}

  //
  // Information about points
  //
  // Get the current point:
  SqrdistIndex<Tdata>* current_point() {
    if (use_landmarks)
      return &point_heap[0];
    else
      return &cur_position[0];
  }
  // Discard the current point and advance to the next-closest:
  void next_point();
  // Test to see if we've exhausted all the points:
  bool is_empty() {
    if (use_landmarks)
      return point_heap_end == point_heap.begin();
    else
      return cur_position == point_heap.end();
  }
};
