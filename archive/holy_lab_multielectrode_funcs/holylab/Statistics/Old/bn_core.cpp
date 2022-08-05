#include <stdexcept>
#include <string.h>
#include <vector>
#include <algorithm>

using namespace std;

// A class to hold the raw data points
template <class dataType>
class data_points {
  public:
    const dataType *x;  // the data points, in order coordinates-then-points
    int d;              // the # of dimensions
    int N;              // the # of points

    data_points() {
      x = NULL;
      d = N = 0;
    }
};

// A class to store information about nearest neighbors
template <class dataType>
class balanced_neighborhood {
  private:
    bool internally_allocated;
  public:
    // "Inputs"
    int d;             // dimensionality
    int N;             // number of points in dataset
    dataType *ybase;   // the "base point"
    dataType *l2;      // the square length scales used to measure distances from the base point (as a default, you can set all to 1)
    // "Outputs"
    int n_valid;        // the # of neighbors that have been set (i.e., placed in order of increasing distance)
    int32_t *neighborLabel; // the index of the points, sorted in order of increasing distance (using a metric scaled by l2) from the base point. Only the first n_valid entries are valid. (This variable is declared as fixed-width to ease Matlab compatibility)
    dataType *neighborDist; // the sorted array of distances: neighborDist[i] = dist(x(:,neighborLabel[i]),y)
    int n_z;        // the # of neighbors required to exceed the balancing criterion by a factor of z (pre-backtracking)
    int n;                  // the # of neighbors in the balanced neighborhood (post-backtracking)
    dataType *nbrhoodMean;     // the center-of-mass of the neighbors in the balanced neighborhood
    dataType *nbrhoodMSDisp;   // the coordinatewise-mean square displacement relative to the base point)
    dataType RMSDist;      // the root-mean-square distance of points in the balanced neighborhood (using the length scales in l2 for the metric)

    balanced_neighborhood() {
      internally_allocated = false;
      ybase = l2 = nbrhoodMean = nbrhoodMSDisp = neighborDist = NULL;
      neighborLabel = NULL;
      n_valid = n_z = n = 0;
    }
    void allocate(int di,int Ni) {
      internally_allocated = true;
      d = di;
      N = Ni;
      if (d>0) {
        ybase = new dataType[d];
        l2 = new dataType[d];
        nbrhoodMean = new dataType[d];
        nbrhoodMSDisp = new dataType[d];
      }
      if (N>0) {
        neighborLabel = new int32_t[N];
        neighborDist = new dataType[N];
      }
    }
    ~balanced_neighborhood() {
      if (internally_allocated) {
        delete[] ybase;
        delete[] l2;
        delete[] nbrhoodMean;
        delete[] nbrhoodMSDisp;
        delete[] neighborLabel;
        delete[] neighborDist;
      }
    }

    void copy(const balanced_neighborhood<dataType> &nbr) {
      memcpy(ybase,nbr.ybase,d*sizeof(dataType));
      memcpy(l2,nbr.l2,d*sizeof(dataType));
      memcpy(nbrhoodMean,nbr.nbrhoodMean,d*sizeof(dataType));
      memcpy(nbrhoodMSDisp,nbr.nbrhoodMSDisp,d*sizeof(dataType));
      memcpy(neighborLabel,nbr.neighborLabel,N*sizeof(int32_t));
      memcpy(neighborDist,nbr.neighborDist,N*sizeof(dataType));
      n_valid = nbr.n_valid;
      n_z = nbr.n_z;
      n = nbr.n;
      RMSDist = nbr.RMSDist;
    }
    void copy_yl_to_base(const dataType *yi,const dataType *l2i) {
      if (ybase != yi)
        memcpy(ybase,yi,d*sizeof(dataType));
      if (l2 != l2i)
        memcpy(l2,l2i,d*sizeof(dataType));
    }
    void copy_yl_to_output(const dataType *yi,const dataType *l2i) {
      memcpy(nbrhoodMean,yi,d*sizeof(dataType));
      memcpy(nbrhoodMSDisp,l2i,d*sizeof(dataType));
    }
    void copy_yl_to_base(const balanced_neighborhood<dataType> &nbr) {
      memcpy(ybase,nbr.ybase,d*sizeof(dataType));
      memcpy(l2,nbr.l2,d*sizeof(dataType));
    }
};

// A structure that holds parameters determining features of the algorithms
struct optionStruct {
  double factor;     // the factor by which the step size should exceed the uncertainty in the step size (in terms of distance, not distance-squared) to be certain the nbrhood is big enough to satisfy the balancing criterion
  int min_to_check;  // the minimum # of points included before testing balancing criterion
  bool backtrack;    // if true, will "backtrack" to point where balancing criterion is first satisfied (with factor = 1)
  int n_min;         // the minimum # of points in the balanced neighborhood
  int iter_max;      // maximum number of iterative updating steps
  bool variable_metric;    // if true, coordinate scales will be updated
  bool protect_zeros;  // if true, variable-metric operations will replace l2 == 0 coordinates with their previous value
  double tol;        // tolerance for determining cycling
  int n_threads;     // not used yet

  optionStruct() {
    factor = 3.0;
    min_to_check = 3;  // d+1 is a better default, but need to know d
    backtrack = true;
    n_min = 2;
    iter_max = 1000;
    variable_metric = true;
    protect_zeros = true;
    tol = 1e-6;
    n_threads = 1;
  }
};

// Temporary storage needed for the MSAMS algorithm
template <class dataType>
class msams_temporary_storage {
public:
 // allocate N+1 of the following (heaps work faster with unit-offset indexing)
  dataType *heapDist;
  int *heapLabel;
  // allocate d of the following
  dataType *dxcum;
  dataType *sumsq;
  dataType *dxcum_saved;
  dataType *sumsq_saved;

  msams_temporary_storage() {
    heapDist = NULL;
    heapLabel = NULL;
    dxcum = NULL;
    sumsq = NULL;
    dxcum_saved = NULL;
    sumsq_saved = NULL;
  }
  void allocate(int d,int N) {
    heapDist = new dataType[N+1];
    heapLabel = new int[N+1];
    dxcum = new dataType[d];
    sumsq = new dataType[d];
    dxcum_saved = new dataType[d];
    sumsq_saved = new dataType[d];
  }
  ~msams_temporary_storage() {
    delete[] heapDist;
    delete[] heapLabel;
    delete[] dxcum;
    delete[] sumsq;
    delete[] dxcum_saved;
    delete[] sumsq_saved;
  }
};

// A class for storing trajectories and testing whether a point has entered a cycle
template <class dataType>
class trajectory_operations {
  private:
    int cycle_start;
    bool internally_allocated;
  public:
    int d;
    dataType *ytraj;
    dataType *l2traj;
    int32_t *ntraj;
    int32_t *nntraj;
    int traj_length;   // the current number of points on the trajectory

    trajectory_operations();
    void allocate(int di,int iter_max);
    ~trajectory_operations();

    void append(const dataType *y,const dataType *l2,int n,int32_t nnbr);
    bool is_cycling(dataType tol);
    void get_maxn(dataType *y,dataType *l2);
/*
    void push_yl(const dataType *y,const dataType *l2);   // must call first
    void push_sn(dataType l2scale,int n);                 // must call second
*/
};

template <class dataType,class intType>
    class l2_operations {
      private:
        vector<int> lgroups;        // holds the key associated with each coordinate
        vector<int> umultiplicity;  // holds the multiplicity of each unique key
        vector<int> ulookup;        // holds a coordinate index associated with each key
        vector<dataType> l2tmp;
        vector<dataType> cumsum;
        vector<dataType> l2min;

        struct greater {
          bool operator()(const dataType v1,const dataType v2) const {
            return v1 > v2;
          }
        };
        dataType maxval(const dataType a,const dataType b) {return (a > b) ? a : b;}
        dataType minval(const dataType a,const dataType b) {return (a < b) ? a : b;}
      public:
        //l2_operations() {;}
        void initialize(const intType *lgroupsi,int di);
        //~l2_operations() {;}
        void set_minimum(const dataType *l2mini);

        void protect_zeros(dataType *l2,const dataType *l2Old);
        void enforce_groups(dataType *l2in);
        dataType calculate_scalefactor(const dataType *l2base,const dataType *l2new);
    };


// Exception classes
class triangle_inequality : public std::runtime_error {
  public:
    triangle_inequality() : std::runtime_error("Error on triangle inequality") { }
};


template <class dataType,class pointIterator<dataType> >
void bn_criterion(pointIterator<dataType> &pI, bn_tmp<dataType> &tmp,dataType z,int &w_z,int &n)
{
  // A "pointIterator" has to provide the following methods:
  //   pI.N(): the total number of points
  //   pI.d(): the # of dimensions
  //   pI.start(): start with the first point
  //   pI++: advance to the next point
  //   pI[dimIndex]: for the current point, the displacement from the base point along dimension dimIndex
  // Note this implements the "rough" criterion
  int dimIterator;
  dataType dx,tmp1,tmp2;
  dataType z2 = z*z;
  tmp.clear();
  for (w_z = 1, pI.start(); w_z <= pI.N(); w_z++,pI++) {
    sq_sum = 0;
    sum_sq = 0;
    for (dimIterator = 0; dimIterator < pI.d(); dimIterator++) {
      dx = pI[dimIterator];
      tmp1 = (tmp.dxsum[dimIterator] += dx);
      tmp2 = (tmp.dx2sum[dimIterator] += dx*dx);
      sq_sum += tmp1*tmp1;
      sum_sq += tmp2;
    }
    if (sq_sum > z2*sum_sq)
      break;
  }
  n = (w_z-(int)z2+1)/2; // +1 for the "ceiling" operation
}

template <class intType>
struct less_index {
  const intType *value;

  less_index(const intType *valueIn) {value = valueIn;}
  void set(const intType *valueIn) {value = valueIn;}

  bool operator()(const int i1,const int i2) const {
    return value[i1] < value[i2];
  }
};

template <class dataType>
class pointsByDistance {
private:
  const dataType *x;   // the data points, in order coordinates-then-points
  int di;              // the # of dimensions
  int Ni;              // the # of points
  dataType *dx;        // displacement from basept in units of ell
  dataType *dxp;       // pointer to current displacement
  bool internally_allocated;
  // the next variables are for using the triangle inequality to triage points
  dataType *yBase;
  dataType *ellBase;
  int32_t *neighborLabelBase;
  dataType *neighborDistBase;
  int *heapLabel;
  int baseIndex;
  int heapIndex;
  less_index<int32_t> li;
public:
  dataType *y;         // location of the probe point
  dataType *ell;       // the length scales used in the metric
  int n_valid;        // the # of neighbors that have been set (i.e., placed in order of increasing distance)
  int32_t *neighborLabel; // the index of the points, sorted in order of increasing distance (using a metric scaled by l2) from the base point. Only the first n_valid entries are valid. (This variable is declared as fixed-width to ease Matlab compatibility)
  dataType *neighborDist; // the array of distances: neighborDist[i] = dist(x(:,i),y). Only the entries cataloged by the first n_valid entries in neighborLabel are valid.

  void startFromScratch();  // Order points without using triangle inequality
  int N() const {return Ni;}
  int d() const {return di;}
  void start();
  void operator++();
  dataType operator[](dimIndex) const {return dxp[dimIndex];}

// This finds the MSAMS-neighborhood of a probe point, and updates the
// probe location and l2 (metric-scaling).
// It is flexible in that it either examines all points, or can use a
// base point and the triangle inequality to limit the number of points
// that need to be examined
template <class dataType,class intType>
int msams_core(const data_points<dataType> &points,
               const balanced_neighborhood<dataType> &nbrbase,
               const dataType *yin,
               const dataType *l2in,
               l2_operations<dataType,intType> &l2func,
               const intType *assignment,
               bool using_base_neighbor_sorting,
               balanced_neighborhood<dataType> &nbr_out,
	       bool &is_assigned,
               const optionStruct &ops,
               msams_temporary_storage<dataType> &tmp)
{
  //Inputs:
  //  yin is the probe point location
  //  l2in is the square-length scaling each coordinate
  //  If you want to use the triangle inequality to avoid having to test all data points, set using_base_neighbor_sorting=true and supply nbrbase (which you can obtain as nbr_out when calling this function once with using_base_neighbor_sorting=false)
  //  assignment is a vector of length N. If zero, it means the corresponding point in x has not been assigned to a target; if greater than zero, this is the index (unit-offset) of the target. Use this to truncate iteration early.
  //  options should be set appropriately
  //  tmp needs to be allocated

  //Outputs:
  // nearestNbrIndex will be set to the index of the nearest neighbor (note this info might also be in nbrList, but if the user set this to NULL then we will need nearestNbrIndex)
  // nmsams will be the number of points needed to satisfy the MSAMS criterion
  // n will be the number of points "in the neighborhood" (will be less than nmsams if backtrack is true)
  // nbrList, if not input as NULL, will contain the index of the <nmsams> neighbors of the probe point, in order of increasing distance
  // nbrDist, if not input as NULL, will contain the distances to the closest <nmsams> neighbors
  // the return value is the # of points "examined" (placed on the heap), which in general will be more than the number used to satisfy the MSAMS criterion (popped off the heap)

  // Declarations needing no initialization
  int dimIterator;
  int thisPointIndex;
  const dataType *yp, *lp;
  dataType *dxcump, *sumsqp;
  dataType rbound, stepsq;
  int heapI, heapII;
  int tmpInt;
  dataType tmpReal;
  dataType dx;
  dataType l2sum;
  dataType lfactor;    // sqrt(factor for converting length scales)
  dataType r0;         // the distance between ybase and y, using l2base
  dataType rcur;       // distance to closest known point
  int n;
  int n_saved;          // used for backtracking
  bool add_more,pop_more;

  // Copy the inputs to nbr_out
  if (yin != NULL && l2in != NULL)
    nbr_out.copy_yl_to_base(yin,l2in);

  if (using_base_neighbor_sorting) {
    // Compute relationship between base point and current point
    // These will be needed for applying the triangle inequality
    r0 = 0;         // the distance between ybase and y, using l2base
    for (dimIterator = 0; dimIterator < points.d; dimIterator++) {
      dx = nbr_out.ybase[dimIterator] - nbrbase.ybase[dimIterator];
      r0 += dx*dx/nbrbase.l2[dimIterator];
    }
    r0 = sqrt(r0);
    dataType S = l2func.calculate_scalefactor(nbrbase.l2,nbr_out.l2);
    lfactor = sqrt(S);
    dataType tol = 1e-6;
    lfactor *= 1+tol;  // to overwhelm roundoff errors in triangle inequality
    r0 *= 1+tol;
    // initialize rcur to guarantee we start by pushing onto heap
    rcur = nbrbase.neighborDist[points.N-1]/lfactor+1;
  }

  // Initialization for the loop
  n = nbr_out.n_z = nbr_out.n_valid = 0;
  int n_points_to_remove = 0;
  const dataType *yend = nbr_out.ybase + points.d;
  const dataType *xp = points.x;
  bool msams_criterion_met = false;
  for (dimIterator = 0; dimIterator < points.d; dimIterator++) {
    tmp.dxcum[dimIterator] = 0;
    tmp.sumsq[dimIterator] = 0;
  }
  dataType sumsqtot = 0;
  dataType factor2 = ops.factor*ops.factor;
  bool backtrack_met = false;
  int baseIndex = 0;
  int heapIndex = 0;
  is_assigned = false;

  //
  // The main loop
  //
  while (!msams_criterion_met && baseIndex < points.N) {
    //
    // Add points to heap.
    // This may be on an as needed-basis to insure (via the triangle
    // inequality) that we know the next-closest point
    if (using_base_neighbor_sorting)
      add_more = (rbound = (nbrbase.neighborDist[baseIndex]-r0)/lfactor) <= rcur;
    else
      add_more = true;
    while (baseIndex < points.N && add_more) {
      // Push the next point on the heap
      //
      // First, compute the distance from y
      if (using_base_neighbor_sorting) {
	thisPointIndex = nbrbase.neighborLabel[baseIndex];
	xp = points.x + points.d*thisPointIndex;
      }
      else
	thisPointIndex = baseIndex;
      tmpReal = 0;   // will hold the distance
      for (yp = nbr_out.ybase, lp = nbr_out.l2; yp < yend; xp++,yp++,lp++) {
        dx = *xp-*yp;
        tmpReal += (dx*dx) / *lp;
      }
      tmpReal = sqrt(tmpReal);
      // Sanity check
      if (using_base_neighbor_sorting && tmpReal < rbound) {
        throw triangle_inequality();
      }
      // Prepare to place on heap
      heapI = heapIndex+1;  // unit offset indexing
      // Re-heapify (propagate up the heap to find its proper place)
      heapII = heapI/2;
      while (heapII > 0 && tmp.heapDist[heapII] > tmpReal) {
        tmp.heapDist[heapI] = tmp.heapDist[heapII];  // shift the bigger point down
        tmp.heapLabel[heapI] = tmp.heapLabel[heapII];
        heapI = heapII;
        heapII = heapII/2;
      }
      tmp.heapDist[heapI] = tmpReal;
      tmp.heapLabel[heapI] = thisPointIndex;
      // Prepare for the next point
      baseIndex++;
      heapIndex++;
      if (using_base_neighbor_sorting) {
	rcur = tmp.heapDist[1];  // this is the base point (the closest non-processed point yet)
	add_more = (rbound = (nbrbase.neighborDist[baseIndex]-r0)/lfactor) <= rcur;
      }
    }


    //
    // The heap is sufficiently loaded that we can read one or more points off
    //
    pop_more = true;
    while (heapIndex > 0 && pop_more) {
      thisPointIndex = tmp.heapLabel[1];  // the closest non-processed point
      if (n == 0) {
        // This is the nearest neighbor
        // Determine whether the nearest neighbor has already been assigned a fate
        if (assignment != NULL && assignment[thisPointIndex])
          is_assigned = true;
      }
      nbr_out.neighborLabel[n] = thisPointIndex;
      nbr_out.neighborDist[n] = tmp.heapDist[1];
      n++;
      if (!msams_criterion_met) {
	// Test MSAMS criterion
	xp = points.x + points.d*thisPointIndex;
	stepsq = 0;
	for (yp = nbr_out.ybase, lp = nbr_out.l2, dxcump = tmp.dxcum, sumsqp = tmp.sumsq; yp < yend; xp++, yp++, lp++, dxcump++, sumsqp++) {
	  dx = *xp-*yp;
	  *dxcump += dx;
	  tmpReal = dx*dx / *lp;
	  sumsqtot += tmpReal;
	  *sumsqp += tmpReal;
	  stepsq += (*dxcump * *dxcump) / *lp;
	}
	msams_criterion_met = n >= ops.min_to_check && stepsq > factor2 * sumsqtot;
        if (msams_criterion_met) {
          nbr_out.n_z = n;
          n_points_to_remove = (int) (stepsq/sumsqtot) + 1;
          if (n_points_to_remove >= n)
            n_points_to_remove = n-ops.n_min;
        }
      }

      // Take point off heap
      tmpReal = tmp.heapDist[heapIndex]; // move last point to top
      tmpInt = tmp.heapLabel[heapIndex];
      heapIndex--;
      heapI = 1; // re-heapify
      heapII = 2;
      dataType tmpReal2,tmpReal3;   // to speed up heap operations
      while (heapII <= heapIndex) {
	tmpReal2 = tmp.heapDist[heapII];
	if (heapII < heapIndex)
	  if (tmpReal2 > (tmpReal3 = tmp.heapDist[heapII+1])) {
	    heapII++;  // choose the smaller of the two children
	    tmpReal2 = tmpReal3;
	  }
	if (tmpReal <= tmpReal2)
	  break; // heap is fine from here on out, quit
	tmp.heapDist[heapI] = tmpReal2;
	tmp.heapLabel[heapI] = tmp.heapLabel[heapII];
	heapI = heapII;
	heapII = 2*heapII;
      }
      tmp.heapDist[heapI] = tmpReal;
      tmp.heapLabel[heapI] = tmpInt;

      // Determine whether we should keep popping points
      if (using_base_neighbor_sorting) {
        // When using neighbors, we only need to satisfy the MSAMS criterion before quitting. But, we may have to add more points to the heap before we can pop any more.
	pop_more = !msams_criterion_met &&
	  (tmp.heapDist[1] <= rbound || baseIndex == points.N);
      }
      else {
	// If we're not using the triangle inequality, we probably
	// want to pop them all. The exception: if the nearest
	// neighbor has been assigned, we can stop after obtaining the MSAMS neighborhood.
	pop_more = !msams_criterion_met || !is_assigned;
      }
    }
    // Prepare for next round of pushing onto heap
    if (using_base_neighbor_sorting)
      if (heapIndex > 0)
	rcur = tmp.heapDist[1];
      else
	rcur = nbrbase.neighborDist[points.N-1]/lfactor+1;  // this will force a push
  }
  //
  // We have satisfied the MSAMS criterion (or we have used all N points)
  //
  nbr_out.n_valid = n;
  // Do backtracking
  if (msams_criterion_met && ops.backtrack) {
    for (n = nbr_out.n_z; n > nbr_out.n_z - n_points_to_remove; n--) {
      thisPointIndex = nbr_out.neighborLabel[n-1];
      xp = points.x + points.d*thisPointIndex;
      yp = nbr_out.ybase;
      lp = nbr_out.l2;
      for (dxcump = tmp.dxcum, sumsqp = tmp.sumsq; yp < yend; xp++, yp++, lp++, dxcump++, sumsqp++) {
        dx = *xp - *yp;
        *dxcump -= dx;
        tmpReal = dx*dx / *lp;
        *sumsqp -= tmpReal;
      }
    }
  }
  // Compute root-mean-square distance
  dataType RMSDist = 0;
  for (dimIterator = 0; dimIterator < points.d; dimIterator++)
    RMSDist += tmp.sumsq[dimIterator];
  nbr_out.RMSDist = sqrt(RMSDist/n);
  // Update the values
  for (dimIterator = 0; dimIterator < points.d; dimIterator++)
    nbr_out.nbrhoodMean[dimIterator] = nbr_out.ybase[dimIterator] + tmp.dxcum[dimIterator]/n;
  for (dimIterator = 0; dimIterator < points.d; dimIterator++)
    nbr_out.nbrhoodMSDisp[dimIterator] = nbr_out.l2[dimIterator] * tmp.sumsq[dimIterator]/n;

  nbr_out.n = n;
  return baseIndex;
}

// This algorithm moves a point repeatedly to the mean of its balanced neighborhood, until it converges to a fixed point or starts cycling
template <class dataType,class intType>
bool msams_converge_core(const data_points<dataType> &points,
                         balanced_neighborhood<dataType> &nbrbase,
                         const dataType *y,
                         const dataType *l2,
                         l2_operations<dataType,intType> &l2func,
                         const intType *assignment,
                         balanced_neighborhood<dataType> &nbr,
                         dataType *minDist,
                         bool &is_assigned,
                         int *niterp,
                         trajectory_operations<dataType> &trajectory,
                         const optionStruct &ops)
{
  int dimIterator,ptIterator;
  int n_examined;
  balanced_neighborhood<dataType> *nbrP;  // points to most-recent nbrInfo
  bool using_base_neighbor_sorting;
  bool updatey,updatel2;

  // Allocation
  msams_temporary_storage<dataType> mtmp;
  mtmp.allocate(points.d,points.N);

  // Initialization
  int n_extra;    // holds the cumulative # of extra points examined since last re-basing
  if (nbrbase.n_valid == points.N)
    n_extra = 0;   // we'll use triangle inequality to try to limit the search
  else
    n_extra = points.N;
  bool isdone = false;
  int n_iter = 0;
  if (niterp != NULL)
    niterp[0] = niterp[1] = 0;
  if (minDist != NULL)
    fill(minDist,minDist+points.N,-1);   // a sentinel to specify hasn't been determined yet
  nbr.copy_yl_to_output(y,l2);  // set the inputs as if they are output from a previous round
  nbrP = &nbr;

  // Do the computation: move 1 point until it stops moving
  // (or starts cycling)
  while (!isdone) {
    // Move the probe point
    l2func.enforce_groups(nbrP->nbrhoodMSDisp);  // (if variable_metric = false, all will be in a single group)
    using_base_neighbor_sorting = (n_extra < points.N/2);
    if (using_base_neighbor_sorting) {
      // Use neighbor information accumulated previously to speed the computation
      n_examined = msams_core<dataType,intType>(points,nbrbase,nbrP->nbrhoodMean,nbrP->nbrhoodMSDisp,l2func,assignment,using_base_neighbor_sorting,nbr,is_assigned,ops,mtmp);
      nbrP = &nbr;
      n_extra += n_examined - nbr.n_z;
      if (niterp)
        niterp[1]++;
    } else {
      // It's time to "re-base" the calculation around the current
      // point (i.e., start from scratch)
      n_examined = msams_core<dataType,intType>(points,nbr,nbrP->nbrhoodMean,nbrP->nbrhoodMSDisp,l2func,assignment,using_base_neighbor_sorting,nbrbase,is_assigned,ops,mtmp);
      nbrP = &nbrbase;
      n_extra = 0;
      if (niterp)
        niterp[0]++;
    }
    if (ops.protect_zeros)
      l2func.protect_zeros(nbrP->nbrhoodMSDisp,nbrP->l2);

    // Cycle-testing
    if (!is_assigned) {
      trajectory.append(nbrP->ybase,nbrP->l2,nbrP->n,nbrP->neighborLabel[0]);
      // Compare this point against previous points to see if we've
      // reached steady-state or a cycle
      isdone = trajectory.is_cycling(ops.tol);
      if (!isdone) {
        n_iter++;
        isdone = n_iter == ops.iter_max;
      }
    }

    // Update minDist
    if (minDist != NULL) {
      dataType rmsd = nbrP->RMSDist;
      if (rmsd > 0) {
        for (ptIterator = 0; ptIterator < nbrP->n_valid; ptIterator++) {
          int thisPt = nbrP->neighborLabel[ptIterator];
          dataType thisDist = nbrP->neighborDist[ptIterator]/rmsd;
          if (minDist[thisPt] < 0 || thisDist < minDist[thisPt])
            minDist[thisPt] = thisDist;
        }
      }
    }

    // If nearest neighbor is already assigned, return
    if (is_assigned) {
      // If our most recent iteration put its output in nbrbase, we have to copy the data to nbr
      if (nbrP == &nbrbase)
        nbr.copy(nbrbase);
      return false;   // didn't converge, but exited because of assignment
    }
  }
  if (n_iter < ops.iter_max) {
    // It entered a cycle. Keep the position with largest n
    trajectory.get_maxn(nbr.ybase,nbr.l2);
    // Since this almost certainly wasn't the most recent point, fill in the appropriate details with a final call to msams_core
    msams_core<dataType,intType>(points,nbrbase,nbr.ybase,nbr.l2,l2func,NULL,true,nbr,is_assigned,ops,mtmp);

    /*
    mexPrintf("ntraj:\n");
    for (dimIterator = 0; dimIterator < n_iter; dimIterator++)
    mexPrintf("  %d\n",ntraj[dimIterator]);
    */
  } else
    nbr.copy(*nbrP);
  return (n_iter < ops.iter_max);
}


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


template <class dataType>
    trajectory_operations<dataType>::trajectory_operations()
{
  ytraj = l2traj = NULL;
  ntraj = nntraj = NULL;
  d = 0;
  traj_length = 0;
  cycle_start = -1;
  internally_allocated = false;
}

template <class dataType>
    void trajectory_operations<dataType>::allocate(int di,int iter_max)
{
  d = di;
  internally_allocated = true;
  ytraj = new dataType[d*iter_max];
  l2traj = new dataType[d*iter_max];
  ntraj = new int32_t[iter_max];
  nntraj = new int32_t[iter_max];
}

template <class dataType>
    trajectory_operations<dataType>::~trajectory_operations()
{
  if (internally_allocated) {
    delete[] ytraj;
    delete[] l2traj;
    delete[] ntraj;
    delete[] nntraj;
  }
}

template <class dataType>
    void trajectory_operations<dataType>::append(const dataType *y,const dataType *l2,int n,int32_t nnbr)
{
  dataType *thisy = ytraj + traj_length*d;
  dataType *thisl = l2traj + traj_length*d;
  for (int dimIterator = 0; dimIterator < d; dimIterator++) {
    thisy[dimIterator] = y[dimIterator];
    thisl[dimIterator] = l2[dimIterator];
  }
  ntraj[traj_length] = n;
  nntraj[traj_length] = nnbr;

  traj_length++;
}

template <class dataType>
    bool trajectory_operations<dataType>::is_cycling(dataType tol)
{
  // Note this function sets cycle_start to the index at which a cycle is first detected (i.e., it lies between cycle_start and the end)
  // Convergence is tested using only "dimensionless" combinations: dy^2/l2 and |diff(l2)|/l2
  const dataType *y = ytraj + (traj_length-1)*d;
  const dataType *l2 = l2traj + (traj_length-1)*d;
  const int cycle_min = 3;   // don't declare convergence unless has been stationary for this long (1 being the minimum)
  cycle_start = traj_length - 1 - cycle_min;
  const dataType *ytrajtmp,*l2trajtmp;
  dataType dy,suml2,dl2,err;
  int dimIterator;
  for (ytrajtmp = ytraj + d*cycle_start,l2trajtmp = l2traj+d*cycle_start; ytrajtmp >= ytraj; ytrajtmp-=d,l2trajtmp-=d,cycle_start--) {
    err = 0;
    for (dimIterator = 0; dimIterator < d; dimIterator++) {
      dy = y[dimIterator]-ytrajtmp[dimIterator];
      suml2 = l2trajtmp[dimIterator]+l2[dimIterator];
      dl2 = l2trajtmp[dimIterator]-l2[dimIterator];
      err += dy*dy/suml2;
      dl2 = (dl2 < 0) ? -dl2 : dl2;  // absolute value
      err += dl2/suml2;
    }
    if (err < tol)
      return true;
  }
  return false;
}

template <class dataType>
    void trajectory_operations<dataType>::get_maxn(dataType *y,dataType *l2)
{
  int n_max = 0;
  int maxIndex = -1;
  for (int cycleIterator = traj_length-1; cycleIterator >= cycle_start; cycleIterator--) {
    if (n_max < ntraj[cycleIterator]) {
      n_max = ntraj[cycleIterator];
      maxIndex = cycleIterator;
    }
  }
  dataType *ytrajtmp, *l2trajtmp;
  int dimIterator;
  for (ytrajtmp = ytraj + d*maxIndex,
       l2trajtmp = l2traj + d*maxIndex,
       dimIterator = 0; dimIterator < d; dimIterator++) {
         y[dimIterator] = ytrajtmp[dimIterator];
         l2[dimIterator] = l2trajtmp[dimIterator];
       }
}

template <class dataType,class intType>
    void l2_operations<dataType,intType>::initialize(const intType *lgroupsi,int di)
{
  // We know that there will be d unique values, so to avoid any need for overflow protection we will re-map to consecutive integers
  lgroups.clear();
  umultiplicity.clear();
  ulookup.clear();
  l2tmp.clear();
  cumsum.clear();
  if (lgroupsi != NULL) {
    // Calculate an index table
    vector<int> sortIndex(di);
    int dimIterator;
    for (dimIterator = 0; dimIterator < di; dimIterator++)
      sortIndex[dimIterator] = dimIterator;
    // rearrange sortIndex so lgroupsi[sortIndex] is sorted
    less_index li(lgroupsi);
    sort(sortIndex.begin(),sortIndex.end(),li);
    // calculate the consecutive-integer key associated with each input key
    int this_tag = 0;
    lgroups.resize(di);
    lgroups[sortIndex[0]] = this_tag;
    ulookup.push_back(sortIndex[0]);
    for (dimIterator = 1; dimIterator < di; dimIterator++) {
      if (li(dimIterator-1,dimIterator)) {
        this_tag++;
        ulookup.push_back(sortIndex[dimIterator]);
      }
      lgroups[sortIndex[dimIterator]] = this_tag;
    }
    // calculate the multiplicity of each unique index
    int n_tags = this_tag+1;
    umultiplicity.resize(n_tags);
    vector<int>::iterator vii;
    for (vii = lgroups.begin(); vii < lgroups.end(); vii++)
      umultiplicity[*vii]++;
    l2tmp.resize(n_tags);
    cumsum.resize(n_tags);
  } else {
    // All the points will be in a single group
    lgroups.resize(di);
    fill(lgroups.begin(),lgroups.end(),0);
    umultiplicity.push_back(di);
    ulookup.push_back(0);
    l2tmp.resize(1);
    cumsum.resize(1);
  }
}

template <class dataType,class intType>
    void l2_operations<dataType,intType>::set_minimum(const dataType *l2mini)
{
  int d = lgroups.size();
  l2min.clear();
  l2min.resize(d);
  for (int i = 0; i < d; i++)
    l2min[i] = l2mini[i];
}

template <class dataType,class intType>
    void l2_operations<dataType,intType>::protect_zeros(dataType *l2,const dataType *l2Old)
{
  int i;
  for (i = 0; i < lgroups.size(); i++)
    if (l2[i] <= 0)
      l2[i] = l2Old[i];
}

template <class dataType,class intType>
    void l2_operations<dataType,intType>::enforce_groups(dataType *l2)
{
  int i;
  fill(l2tmp.begin(),l2tmp.end(),0);
  // Calculate sum in each group
  for (i = 0; i < lgroups.size(); i++)
    l2tmp[lgroups[i]] += l2[i];
  // Calculate the mean in each group
  for (i = 0; i < l2tmp.size(); i++)
    l2tmp[i] /= umultiplicity[i];
  // Set l2 to the mean across the group
  for (i = 0; i < lgroups.size(); i++)
    l2[i] = l2tmp[lgroups[i]];
  // Prevent from being less than the minimum
  for (i = 0; i < l2min.size(); i++)
    if (l2[i] < l2min[i])
      l2[i] = l2min[i];
}

template <class dataType,class intType>
    dataType l2_operations<dataType,intType>::calculate_scalefactor(const dataType *l2base,const dataType *l2new)
{
  int i;
  int n = l2tmp.size();
  for (i = 0; i < n; i++)
    l2tmp[i] = l2new[ulookup[i]]/l2base[ulookup[i]];
  sort(l2tmp.begin(),l2tmp.end(),greater());
  cumsum[0] = umultiplicity[0]*l2tmp[0];
  for (i = 1; i < n; i++)
    cumsum[i] = cumsum[i-1]+umultiplicity[i]*l2tmp[i];
  dataType S;
  S = maxval(l2tmp[0],cumsum[n-1]-cumsum[0]);
  for (i = 1; i < n; i++)
    S = minval(S,maxval((1-i)*l2tmp[i] + cumsum[i-1],cumsum[n-1] - cumsum[i]));
  return S;
}
