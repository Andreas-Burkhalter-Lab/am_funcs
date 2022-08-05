#include <unistd.h>  // for # of CPUs

// A structure that is used to control the behavior of the algorithm
struct optionStruct {
  int n_min;
  double factor;
  char terminate_mode;
  bool backtrack;
  bool any_coordinate;
  double convergence_thresh;
  int max_iter;
  int n_threads;
  int n_cpus;

  optionStruct() {
    n_min = 3;
    factor = 3;
    terminate_mode = 'c';
    backtrack = true;
    any_coordinate = false;
    convergence_thresh = 1e-12;
    max_iter = 1000;
    n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    n_threads = n_cpus;
  }
};

// An output structure
template <class Tdata,class Tint>
struct outputStruct {
  Tdata *y;
  Tint *closestDataIndex;
  Tint *n;
  Tdata *R2;
  Tint *n_iter;
  Tint *convergedFlag;
  Tint *ntraj;
  Tdata *ytraj;
  Tint *xnbrI;

  outputStruct() {
    y = NULL;
    closestDataIndex = NULL;
    n = NULL;
    R2 = NULL;
    n_iter = NULL;
    convergedFlag = NULL;
    ntraj = NULL;
    ytraj = NULL;
    xnbrI = NULL;
  }
};

// A catch-all structure that is needed to support multithreading
template <class Tdata,class Tint>
struct thread_data_type {
  const Tdata *x;
  int d;
  int N;
  landmarkStruct<Tdata,Tint> *lminfo;
  int q;
  outputStruct<Tdata,Tint> *out;
  const optionStruct *ops;

  int probePtIndex_start;
  int probePtIndex_end;
};

/*
// A comparison function for heapsort routines, arranged so that the
// top of the heap will be the closest point
struct is_farther: public binary_function<const SqrdistIndex&,const SqrdistIndex&,bool> {
  bool operator()(const SqrdistIndex &sd1,const SqrdistIndex &ds2) {return sd1.R2 > sd2.R2;}};
*/

template <class Tdata,class Tint>
extern int msams_core(const Tdata *x,int d,int N,landmarkStruct<Tdata,Tint> &lminfo,const Tdata *y,int q,outputStruct<Tdata,Tint> &out,const optionStruct &ops);
template <class Tdata,class Tint>
extern void msams_core_work(thread_data_type<Tdata,Tint> *td);


template <class Tdata,class Tint>
extern void expand_neighborhood_msams(const Tdata *x,int d,const Tdata *thisY,landmarked_neighbors<Tdata,Tint> &lm_nbrs,const optionStruct *ops,vector<Tdata> &x_cum,vector<Tdata> &deltax_cum,vector<Tdata> &sumsq_cum,vector<Tdata> &x_backtrack,int &n,Tdata &R2,int offset,Tint *pointIndex);
