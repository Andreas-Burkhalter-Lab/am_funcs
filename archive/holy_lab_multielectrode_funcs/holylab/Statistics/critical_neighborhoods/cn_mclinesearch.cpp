#ifdef MAIN
#include "mat.h"
#include "MatlabIO.h"
#include "MatlabStubs.h"
#else
#include "mex.h"
#endif

#include "MonteCarloLineSearch.cpp"
#include "BoostThreadPool.h"

// Copyright 2011 by Timothy E. Holy

// Syntax:
//   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim)
//   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim,npoints)
//   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim,npoints,nthreads)

// The following bits of code are needed for multitasking

// The task structure
struct Task {
  double *pT2;
  double *pT2ls;
  int nsim;
};

// This class defines the function that will run in each thread.
template <typename TaskListType,typename T2Type>
class ThreadHandler {
public:
  // The constructor
  // Each thread needs to use a different random seed, or otherwise
  // all threads will use the same random sequence
  ThreadHandler(TaskListType &tasklist,const std::vector<int>& nlist,int d,int npoints,T2Type& T2D,int seed) : mrTaskList(tasklist), mrNList(nlist), mD(d), mNPoints(npoints), mT2D(T2D) { mRng.seed(seed); }

  // The function running in each thread
  void operator()() {
    typename TaskListType::TaskType T;
    bool havetask;

    while (true) {
      // Check if there is another task, and get it. Note these have
      // to be unified into a single call to avoid race conditions,
      // hence the syntax of the "next" member function.
      havetask = mrTaskList.next(T);
      if (!havetask)
        break;  // no more tasks remain
      // Execute the task
      MonteCarloLineSearch(T.pT2,T.pT2ls,T.nsim,mrNList,mD,mNPoints,mT2D,mRng);
    }
  }

private:
  TaskListType&     mrTaskList;
  const std::vector<int>& mrNList;
  int               mD;
  int               mNPoints;
  T2Type            mT2D;
  boost::mt19937    mRng;
};



template <typename T2Type>
void mexThreadHandler(double *pT2,double *pT2ls,int nsim,const std::vector<int>& nlist,int d,int npoints,T2Type &T2D,int nthreads)
{
  if (nthreads == 0) {
    // Run without threads (also not interruptable, no progress report)
    boost::mt19937 rng;
    MonteCarloLineSearch(pT2,pT2ls,nsim,nlist,d,npoints,T2D,rng);
  } else {
    // Multi-task, and perhaps multithreaded.
    // Make it display progress & be interruptable
    typedef TaskList<Task> TaskListType;
    TaskListType tasks(true,true);
    // Split the simulations up into tasks
    int i;
    int ntasks = 100;  // the approximate # of tasks we desire, may be different
    int nsim_per_task = nsim/ntasks;
    if (nsim_per_task < 1) {
      nsim_per_task = nsim/nthreads;
      if (nsim_per_task < 1)
	nsim_per_task = 1;
    }
    Task T;
    T.pT2 = pT2;
    T.pT2ls = pT2ls;
    T.nsim = nsim_per_task;
    double* pT2end = pT2 + nsim*nlist.size();
    int skip_per_task = nsim_per_task*nlist.size();
    while (T.pT2 < pT2end) {
      tasks.push(T);  // add task to the queue
      //mexPrintf("pT2 %d, nsim %d\n",T.pT2,T.nsim);
      T.pT2 += skip_per_task;
      T.pT2ls += skip_per_task;
      if (T.pT2+skip_per_task > pT2end) {
	nsim_per_task = (pT2end-T.pT2)/nlist.size();
	skip_per_task = nsim_per_task*nlist.size();
      }
      T.nsim = nsim_per_task;
    }
    //mexPrintf("%d tasks\n",tasks.size());

    typedef ThreadHandler<TaskListType,T2Type> ThreadHandlerType;
    if (nthreads == 1) {
      // Single-threaded, run in main thread
      // This tests the task-based part without launching new threads
      ThreadHandlerType single(tasks,nlist,d,npoints,T2D,0);
      single();  // This runs all the tasks in the main thread
    } else {
      // Multithreaded
      // Launch the threads. Note we also let the main thread handle a
      // set of tasks---this is necessary if you want to display progress,
      // because that can be done safely only from the main thread.
      boost::thread_group threads;
      for (i = 1; i < nthreads; i++) {  //starting i=1 is deliberate, see above
	ThreadHandlerType multi(tasks,nlist,d,npoints,T2D,i);
	threads.create_thread(multi);
      }
      ThreadHandlerType multi(tasks,nlist,d,npoints,T2D,0); // Set up the main thread processing
      multi();                          // Get the main thread going
      // Wait for all threads to finish
      threads.join_all();
    }
  }
}


void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  double *p,*pT2,*pT2ls;
  int d,nsim,npoints,nthreads,i;
  std::vector<int> nlist;

  if (nrhs < 4 || nrhs > 6)
    mexErrMsgTxt("cn_mclinesearch: requires four to six inputs");
  if (nlhs != 2)
    mexErrMsgTxt("cn_mclinesearch: requires two outputs");

  // Parse the input
  // d
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("cn_mclinesearch: d must be a real scalar");
  d = (int) mxGetScalar(curarg);

  // nlist
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("cn_mclinesearch: nlist must be a real vector");
  p = mxGetPr(curarg);
  for (i = 0; i < mxGetNumberOfElements(curarg); i++)
    nlist.push_back((int) p[i]);

  // covarianceModel
  int covariance_model = 0;
  curarg = prhs[2];
  if (!mxIsChar(curarg))
    mexErrMsgTxt("covarianceModel must be a string");
  mxChar *pc = mxGetChars(curarg);
  if (*pc == 'i')
    covariance_model = 0;
  else if (*pc == 'd')
    covariance_model = 1;
  else if (*pc == 'f')
    covariance_model = 2;
  else
    mexErrMsgIdAndTxt("cn:covariance","covarianceModel %s not recognized",pc);

  // nsim
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("cn_mclinesearch: nsim must be a real scalar");
  nsim = (int) mxGetScalar(curarg);

  // npoints
  int nlist_max = *std::max_element(nlist.begin(),nlist.end());
  npoints = nlist_max * 2 * d;
  if (nrhs > 4) {
    curarg = prhs[4];
    if (!mxIsEmpty(curarg)) {
      if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
	mexErrMsgTxt("cn_mclinesearch: npoints must be a real scalar");
      npoints = (int) mxGetScalar(curarg);
    }
  }
  
  // nthreads
  // Determine the default number of threads based on the number of cores
  nthreads = boost::thread::hardware_concurrency();
  if (nthreads > 2)
    nthreads--;  // save one core (to be nice)
  if (nthreads > 8)
    nthreads = 8;
  if (nrhs > 5) {
    curarg = prhs[5];
    if (!mxIsEmpty(curarg)) {
      if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
	mexErrMsgTxt("cn_mclinesearch: nthreads must be a real scalar");
      nthreads = (int) mxGetScalar(curarg);
    }
  }
  
  mexPrintf("d = %d, length(nlist) = %d, covModel = %d, nsim = %d, npoints = %d, nthreads = %d\n",d,nlist.size(),covariance_model,nsim,npoints,nthreads);

  // Allocate output storage
  plhs[0] = mxCreateDoubleMatrix(nlist.size(),nsim,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nlist.size(),nsim,mxREAL);
  pT2 = mxGetPr(plhs[0]);
  pT2ls = mxGetPr(plhs[1]);

  // Run the algorithm
  if (covariance_model == 0) {
    T2Direct<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>,Moments::Isotropic> T2D(d,nlist_max);
    mexThreadHandler(pT2,pT2ls,nsim,nlist,d,npoints,T2D,nthreads);
  } else if (covariance_model == 1) {
    T2Direct<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>,Moments::Diagonal> T2D(d,nlist_max);
    mexThreadHandler(pT2,pT2ls,nsim,nlist,d,npoints,T2D,nthreads);
  } else {  
    T2Direct<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>,Moments::Full> T2D(d,nlist_max);
    mexThreadHandler(pT2,pT2ls,nsim,nlist,d,npoints,T2D,nthreads);
  }
}

#ifdef MAIN
int main(int argc,char *argv[]) {
//   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim)
  const mxArray *argsIn[4];
  MatlabIO::Load reader(argv[1]);
  argsIn[0] = reader.load("d");
  argsIn[1] = reader.load("nlist");
  argsIn[2] = reader.load("covarianceModel");
  argsIn[3] = reader.load("nsim");

  mxArray *argsOut[2];
  
  mexFunction(2,argsOut,4,argsIn);

  mxDestroyArray(argsOut[0]);
  mxDestroyArray(argsOut[1]);
}
#endif // MAIN
