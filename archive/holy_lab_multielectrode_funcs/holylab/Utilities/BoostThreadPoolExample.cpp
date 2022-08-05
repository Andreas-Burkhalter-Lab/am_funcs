// Here's a complete simple example showing how to use
// BoostThreadPool. Each "task" is simply to sleep for the specified
// amount of time. For a more sophisticated example, see
// BoostThreadPoolExample2.cpp.

// Matlab MEX version: compile with
//    mexpp -boostthread BoostThreadPoolExample.cpp 
// for dialog-based interrupts or
//    mexpp -boostthread -ctrlc BoostThreadPoolExample.cpp 
// if you want to use CTRL-C to interrupt execution.
// For a MAT standalone application, use
//    mexpp -standalone -boostthread BoostThreadPoolExample.cpp

// For a version entirely independent of  Matlab, compile with
//   g++ -DNOMAT BoostThreadPoolExample.cpp -o BoostThreadPoolExample -lboost_thread

#if defined(MAIN)
  // Build a standalone MAT function (with mexpp -standalone ...)
  #include "mat.h"
  #include "MatlabStubs.h"
#elif defined(NOMAT)
  // Build a completely standalone function (with g++ ...)
#else
  // Build a MEX file (with mexpp ...)
  #include "mex.h"
#endif

#include <sstream>
#include "BoostThreadPool.h"


// This task structure holds whatever task-specific information is
// required (in general, it might hold task-specific inputs, a pointer
// to where to store the output, etc.)
struct Task {
  double    sleepTimeInSeconds;
};


// This class defines the function that will run in each thread; when
// you create an object obj of this type, then obj() launches
// processing.
template <typename TaskListType>
class ThreadHandler {
public:
  // The constructor
  ThreadHandler(TaskListType &tasklist) : mrTaskList(tasklist) {;}

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
      //
      // Execute the task.
      //
      // For an existing function that you want to adapt to a
      // task-based workflow, at this point you could call your
      // single-threaded "work function."  However, for new projects
      // it's probably better to have the actual computation inside
      // this while loop.  The main advantage is the ability to call
      // mrTaskList.queueMessage(), so that you can display messages
      // (e.g., for debugging).
      // Note that mex* and mat* functions (and often stdio/iostream)
      // are generally not thread-safe, and should be avoided at all
      // costs.
      // "Display" a message
      {
	std::ostringstream ss;
	ss << "Thread " << boost::this_thread::get_id() << " is sleeping for " << T.sleepTimeInSeconds << " second(s).";
	mrTaskList.queueMessage(ss.str());
      }
      // Execute the task
      boost::this_thread::sleep(boost::posix_time::seconds(T.sleepTimeInSeconds));
    }
  }  // end of task-execution function
private:
  TaskListType&     mrTaskList;
};


// Here's the function that breaks the whole problem up into
// tasks. Since it appears that this will be highly problem-specific,
// no generic solution to the problem is provided; you have to do this
// part yourself. Once the tasks are created, launching the threads
// is easy.
void mexWorkFunctionThreaded(double totalSleepTimeInSeconds,int ntasks,int nthreads)
{
  // Make it display progress & be interruptable
  typedef TaskList<Task> TaskListType;
  TaskListType tasks(true,true);
  // Break the problem into tasks & queue them
  // (this problem is unusually straightforward to break up because
  // there is no output)
  int i;
  Task T;
  T.sleepTimeInSeconds = totalSleepTimeInSeconds/ntasks;
  for (i = 0; i < ntasks; i++)
    tasks.push(T);
  // done breaking the problem into tasks

  // Now launch processing
  typedef ThreadHandler<TaskListType> ThreadHandlerType;
  if (nthreads == 1) {
    // Single-threaded, run in main thread
    // This tests the task-based part without launching new threads
    ThreadHandlerType single(tasks);
    single();  // This runs all the tasks in the main thread
  } else {
    // Multithreaded
    // Launch the threads. Note we also let the main thread handle a
    // set of tasks---this is necessary if you want to display progress,
    // because that can be done safely only from the main thread.
    boost::thread_group threads;
    for (i = 1; i < nthreads; i++) {  //starting i=1 is deliberate, see above
      ThreadHandlerType multi(tasks);
      threads.create_thread(multi);
    }
    ThreadHandlerType multi(tasks); // Set up the main thread processing
    multi();                        // Get the main thread going
    // Wait for all threads to finish
    threads.join_all();
  }
}


// Finally, mexFunction
// The syntax is:
//    BoostThreadPoolExample(totalSleepTimeInSeconds,ntasks)
// or
//    BoostThreadPoolExample(totalSleepTimeInSeconds,ntasks,nthreads)
// if you want to explicitly control the number of threads.
//
// For example:
//    BoostThreadPoolExample(10,10);
// or
//    BoostThreadPoolExample(10,10,1);
// for single-threaded execution.
#ifndef NOMAT
void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  int ntasks,nthreads;
  double totalSleepTimeInSeconds;

  if (nrhs < 2 || nrhs > 3)
    mexErrMsgTxt("BoostThreadPoolExample: requires two or three inputs");

  // Parse the input arguments
  // totalSleepTimeInSeconds
  curarg = prhs[0];
  if (!mxIsDouble(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("BoostThreadPoolExample: totalSleepTimeInSeconds must be a real double scalar");
  totalSleepTimeInSeconds = mxGetScalar(curarg);

  // ntasks
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("BoostThreadPoolExample: ntasks must be an integer");
  ntasks = (int) mxGetScalar(curarg);

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
	mexErrMsgTxt("BoostThreadPoolExample: nthreads must be a real scalar");
      nthreads = (int) mxGetScalar(curarg);
    }
  }
  
  mexPrintf("totalSleepTimeInSeconds %g, ntasks %d, nthreads %d\n",totalSleepTimeInSeconds,ntasks,nthreads);

  // Run the operation
  mexWorkFunctionThreaded(totalSleepTimeInSeconds,ntasks,nthreads);
}  // mexFunction
#endif // ifndef NOMAT




//
// Standalone function
//
// This is helpful for debugging.
#if defined(MAIN) || defined(NOMAT)
int main() {
  double totalSleepTimeInSeconds = 10.0;
  int ntasks = 10;
  int nthreads = 2;

  // You can avoid the following #ifdef by including MatlabStubs.h
  // even in the NOMAT version, but for illustration purposes...
#ifdef MAIN
  mexPrintf("totalSleepTimeInSeconds %g, ntasks %d, nthreads %d\n",totalSleepTimeInSeconds,ntasks,nthreads);
#else
  std::cout << "totalSleepTimeInSeconds " << totalSleepTimeInSeconds << ", ntasks " << ntasks << ", nthreads " << nthreads << std::endl;
#endif

  mexWorkFunctionThreaded(totalSleepTimeInSeconds,ntasks,nthreads);
}
#endif  // defined(MAIN) || defined(NOMAT)
