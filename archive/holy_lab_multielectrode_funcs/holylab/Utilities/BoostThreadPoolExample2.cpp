// Here's a complete simple example showing how to use
// BoostThreadPool.  The elementary calculation is to sum a block of
// numbers (of length "blocksize") and add a random number.  A "task"
// is a set of such blocks, with one output generated for each
// "blocksize" inputs.  At any given time, each thread will be working
// on one task; when the task finishes, the thread checks to see if
// there are any more tasks remaining. If not (or if the user presses
// CTRL-C), the thread finishes.

// Compile with:
//    mexpp -boostthread BoostThreadPoolExample2.cpp
// for dialog-based interrupts or
//    mexpp -boostthread -ctrlc BoostThreadPoolExample2.cpp 
// if you want to use CTRL-C to interrupt execution.

#ifdef MAIN
// Build a standalone function (use mexpp -standalone ...)
#include "mat.h"
#include "MatlabStubs.h"
#else
// Build a standard MEX file
#include "mex.h"
#endif

#include "BoostThreadPool.h"
#include "boost/random.hpp"
#include "boost/random/normal_distribution.hpp"

// First, let's write the single-threaded code, i.e., the code that
// implements your core algorithm. You'd have to write this anyway.
// This part has nothing to do with BoostThreadPool.

// Define a random number generator that lets you get the next
// random number by calling "rng()" and set the seed with
// "rng.seed(int)".
class MyRandomGenerator {
public:
  MyRandomGenerator() : mNormalDistribution(0.0,1.0), mRandomNormal(mRng,mNormalDistribution) {;}
    
  double operator()() {return mRandomNormal();}
  void seed(int i) {mRng.seed(i);}
private:
  boost::mt19937 mRng;
  boost::normal_distribution<> mNormalDistribution;
  boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > mRandomNormal;
};

// This is the "single-threaded work function," which sums numbers in
// blocks and adds a random number.
template <typename RandomGeneratorType>
void mexWorkFunction(double *output,const double *input,int nblocks,int blocksize,RandomGeneratorType &rng)
{
  double tmp;
  int i;
  const double *inputend;

  // Loop over blocks
  for (i = 0; i < nblocks; i++,output++) {
    // First, sum "blocksize" numbers
    tmp = 0;
    inputend = input+blocksize;
    for (; input < inputend; input++)
      tmp += *input;
    // Now add a random number
    tmp += rng();
    // Store the output
    *output = tmp;
  }
}

// Below this point, all the remaining code centers around
// multithreading and task-handling.

// This task structure documents where to store the result, where to
// get the data, and how many blocks there are in this task.
struct Task {
  double       *output;
  const double *input;
  int          nblocks;
};

// This class defines the function that will run in each thread.
template <typename TaskListType,typename RandomGeneratorType>
class ThreadHandler {
public:
  // The constructor
  // Each thread needs to use a different random seed, or otherwise
  // all threads will use the same random sequence
  ThreadHandler(TaskListType &tasklist,int blocksize,int seed) : mrTaskList(tasklist), mBlockSize(blocksize) { mRng.seed(seed); }

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
      mexWorkFunction(T.output,T.input,T.nblocks,mBlockSize,mRng);
      // Pause for a short while, to make sure that we don't have one
      // thread do all the work before others have even started
      // (you would delete this from any production code)
      boost::this_thread::sleep(boost::posix_time::millisec(50));
    }
  }

private:
  TaskListType&     mrTaskList;
  int               mBlockSize;
  MyRandomGenerator mRng;
};

// Here's the function that breaks the whole problem up into
// tasks. Since it appears that this will be highly problem-specific,
// no generic solution to the problem is provided, you have to do this
// part yourself.
void mexWorkFunctionThreaded(double *output,const double *input,int inputlength,int blocksize,int nthreads)
{
  if (nthreads == 0) {
    // Here we've chosen to make this a special case for
    // single-threaded execution that bypasses all the fancy
    // stuff. This can be a good idea as a work-around for cases where
    // the user is getting crashes with the multithreaded version.
    // This won't display progress, nor will it be interruptable.
    MyRandomGenerator rng;
    rng.seed(0);
    mexWorkFunction(output,input,inputlength/blocksize,blocksize,rng);
  } else {
    // Multi-task, and perhaps multithreaded.
    // Make it display progress & be interruptable
    typedef TaskList<Task> TaskListType;
    TaskListType tasks(true,true);
    // Break the problem into tasks
    int i;
    int ntasks = 100;  // the approximate # of tasks we desire, may be different
    int nblocks = inputlength/(ntasks*blocksize); // the size of 1 task
    if (nblocks < 1)
      nblocks = 1;
    Task T;  // prepare the task data
    T.output = output;
    T.input = input;
    T.nblocks = nblocks;
    const double *inputend = input+inputlength;
    while (T.input < inputend) {
      //mexPrintf("input %x, output %x, nblocks %d\n",T.input,T.output,T.nblocks);
      tasks.push(T);  // add a task to the queue
      T.output += nblocks;
      T.input += nblocks*blocksize;
      if (T.input+blocksize*nblocks > inputend)
	nblocks = (inputend-T.input)/blocksize;  // a "fractional-sized" task at the end
      T.nblocks = nblocks;
    }
    // done splitting it up into tasks (whew!)

    typedef ThreadHandler<TaskListType,MyRandomGenerator> ThreadHandlerType;
    if (nthreads == 1) {
      // Single-threaded, run in main thread
      // This tests the task-based part without launching new threads
      ThreadHandlerType single(tasks,blocksize,0);
      single();  // This runs all the tasks in the main thread
    } else {
      // Multithreaded
      // Launch the threads. Note we also let the main thread handle a
      // set of tasks---this is necessary if you want to display progress,
      // because that can be done safely only from the main thread.
      boost::thread_group threads;
      for (i = 1; i < nthreads; i++) {  //starting i=1 is deliberate, see above
	ThreadHandlerType multi(tasks,blocksize,i);
	threads.create_thread(multi);
      }
      ThreadHandlerType multi(tasks,blocksize,0); // Set up the main thread processing
      multi();                          // Get the main thread going
      // Wait for all threads to finish
      threads.join_all();
    }
  }
}


// Finally, mexFunction
// The syntax is:
//    output = BoostThreadPoolExample2(input,blocksize)
// or
//    output = BoostThreadPoolExample2(input,blocksize,nthreads)
// if you want to explicitly control the number of threads.
//
// For example:
//    output = BoostThreadPoolExample2(zeros(1,50),5);
void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  double *pInput, *pOutput;
  int inputlength,blocksize,nthreads;

  if (nrhs < 2 || nrhs > 3)
    mexErrMsgTxt("BoostThreadPoolExample2: requires two or three inputs");

  // Parse the input arguments
  // "input"
  curarg = prhs[0];
  if (!mxIsDouble(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("BoostThreadPoolExample2: input must be real double");
  pInput = mxGetPr(curarg);
  inputlength = mxGetNumberOfElements(curarg);

  // blocksize
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("BoostThreadPoolExample2: blocksize must be an integer");
  blocksize = (int) mxGetScalar(curarg);
  int outputlength = inputlength/blocksize;
  if (outputlength*blocksize != inputlength)
    mexErrMsgTxt("BoostThreadPoolExample2: the length of the input is not a multiple of the blocksize");

  // nthreads
  // Determine the default number of threads based on the number of cores
  nthreads = boost::thread::hardware_concurrency();
  if (nthreads > 2)
    nthreads--;    // be nice and save one core for other uses
  if (nthreads > 8)
    nthreads = 8;  // limit to 8 threads max
  if (nrhs > 2) {
    curarg = prhs[2];
    if (!mxIsEmpty(curarg)) {
      if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
	mexErrMsgTxt("BoostThreadPoolExample2: nthreads must be a real scalar");
      nthreads = (int) mxGetScalar(curarg);
    }
  }
  
  mexPrintf("inputlength %d, outputlength %d, blocksize %d, nthreads %d\n",inputlength,outputlength,blocksize,nthreads);

  // Allocate output storage
  plhs[0] = mxCreateDoubleMatrix(1,outputlength,mxREAL);
  pOutput = mxGetPr(plhs[0]);

  // Run the operation
  mexWorkFunctionThreaded(pOutput,pInput,inputlength,blocksize,nthreads);
}

#ifdef MAIN
int main() {
  const int blocksize = 5;
  const int inputlength = 50;
  const int outputlength = inputlength/blocksize;
  double *pOutput = new double[outputlength];
  double *pInput = new double[inputlength];
  int nthreads = 2;

  mexPrintf("inputlength %d, outputlength %d, blocksize %d\n",inputlength,outputlength,blocksize);

  mexWorkFunctionThreaded(pOutput,pInput,inputlength,blocksize,nthreads);

  delete[] pOutput;
  delete[] pInput;
}
#endif  // MAIN
