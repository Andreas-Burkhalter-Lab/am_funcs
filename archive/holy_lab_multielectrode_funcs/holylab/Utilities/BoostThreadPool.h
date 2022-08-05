#ifndef BOOSTTHREADPOOL_H
#define BOOSTTHREADPOOL_H

//
// Using this class:
//   Constructor:
//     TaskList<TaskType> tasks(bool isShowProgress,bool isInterruptable)
//   If isShowProgress is true, command-line output will be generated
//   indicating how many tasks remain to be processed.  If
//   isInterruptable is true, then the MEX file execution can be
//   interrupted, either by CTRL-C or by dismissing a dialog box.
//
//   Operations:
//     void push(const TaskType &T)
//   Adds a task T to the queue.  Typically, all tasks are added
//   before any threads are launched to process them, but the call is
//   thread-safe so this is not required.
//
//     bool next(TaskType &T)
//   Pop a task from the queue. The new task is deposited in T. The
//   boolean output is true if there is a task that needs processing,
//   and false otherwise. This call is thread-safe.
//
//     void queueMessage(const std::string &s)
//   A thread-safe replacement for iostream/stdio/mexPrintf output.
//   Any accumulated messages are displayed each time the main thread
//   (the one that constructed the TaskList) calls "next", or when the
//   TaskList destructor is called.
//
//     int progressInterval()
//     progressInterval(int i)
//   Get and set, respectively, the minimum number of tasks completed
//   between progress reports. Default value of 1. Increase this to
//   reduce the amount of command-line output when the number of tasks
//   is large.
//
//     int size()
//   Returns the total number of tasks on the queue. This call is not
//   thread-safe, in the sense that the size can change between the
//   time of the call and the time you do something with it.
//
//     bool isMainThread()
//   True if this thread is the one that created the TaskList. This
//   can be used to make thread-unsafe calls specifically from the
//   main thread.
//
// You need to write a certain amount of supporting code, for example
// the code that splits the complete job up into tasks.  Usually this
// step is the main effort of writing mulithreaded code.
//
// See BoostThreadedPoolExample.cpp for a complete worked example.


// When running within Matlab: Boost.Thread should be linked
// statically, see explanation in mexpp.m. If you want CTRL-C
// interrupt support, compile with -DCTRLC -lut. Compiling via mexpp
// takes care of these issues for you.



#include <queue>
#include <string>
#ifndef mex_h
#include <iostream>
#endif
#include "boost/thread.hpp"

#if defined(mex_h) && defined(CTRLC)
// Handling CTRL-C within Malab
// We must declare the function for checking for CTRL-C because it's
// "undocumented"
// Thanks to Wotao Yin, http://www.caam.rice.edu/~wy1/links/mex_ctrl_c_trick/
#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif // __cplusplus
#endif // mex_h && CTRLC

// The class that implements the task list. This is the main work
// needed to support multithreading, and also enables
// progress/interrupt support even when running single-threaded.
template <typename taskType>
class TaskList {
public:
  // Typedefs
  typedef taskType TaskType;

  // The constructor
  TaskList(bool isShowProgress,bool isInterruptable) : mIsShowProgress(isShowProgress), mIsInterruptable(isInterruptable), mIsRunning(true) {
    // Store the thread ID of the creator so that we can ensure that
    // we only display output within the main thread.
    mMainThreadID = boost::this_thread::get_id();
    // Set up progress reporting
    mProgressInterval = 1;  // default minimum # tasks between progress reports
    mLastReport = 1; // initialize to 1 to ensure that the first time through generates a report (unless mProgressInterval gets set higher)
    if (isShowProgress) {
#ifdef mex_h
      mexPrintf("Tasks remaining: ");
#else
      std::cout << "Tasks remaining: " << std::flush;
#endif
    }
    // Set up "interrupts"
#if defined(mex_h) && !defined(CTRLC)
    if (isInterruptable) {
      // Within Matlab, an alternative to libut/CTRL-C for
      // interrupting execution: create a dialog that can be checked
      // for dismissal, signaling the desire to terminate.
      const mxArray *retTrap;
      mxArray *str = mxCreateString(mexFunctionName());
      retTrap = mexCallMATLABWithTrap(1,&mInterruptHandle,1,&str,"mex_interrupt");
      if (retTrap != NULL)
	mexErrMsgTxt("Could not create dialog box for interrupt, aborting");
      mxDestroyArray(str);
    }
#endif
  }
  // The destructor
  ~TaskList() {
    // Display all remaining messages
    bool isMessageShown = mIsShowProgress || !mMessageQueue.empty();
    while (!mMessageQueue.empty()) {
#ifdef mex_h
      mexPrintf("%s\n",mMessageQueue.front().c_str());
#else
      std::cout << mMessageQueue.front() << std::endl;
#endif
      mMessageQueue.pop();
    }
    if (mIsShowProgress)
      if (mTaskQueue.empty()) {
#ifdef mex_h
	mexPrintf("done.\n"); 
#ifndef CTRLC
	// Close the dialog box we opened
	if (mIsInterruptable) {
	  const mxArray *retTrap,*ret;
	  retTrap = mexCallMATLABWithTrap(0,NULL,1,&mInterruptHandle,"delete");
	  if (retTrap != NULL)
	    mexErrMsgTxt("Could not delete dialog box for interrupt, aborting");
	}
#endif // CTRLC
#else
	std::cout << "done." << std::endl;
#endif // mex_h
      }
      else { // mTaskQueue.empty()
#ifdef mex_h	
	mexPrintf("interrupted.\n");
#else
	std::cout << "interrupted." << std::endl;
#endif
      }
#ifdef mex_h
    if (isMessageShown)
      mexEvalString("drawnow;");  // flush screen output
#endif
  }

  // Operations
  void push(const TaskType& T) {
    // Under normal usage patterns, we don't have to worry about
    // locking during "push", because push will be called from the
    // main thread before any of the worker threads are created.
    // Nevertheless, it seems wise to guard against "abuse."
    mTaskMutex.lock();
    mTaskQueue.push(T);
    mLastReport++;
    mTaskMutex.unlock();
  }
  bool next(TaskType& T) {
    // The main principle here is to only hold locks when necessary
#ifdef mex_h
    // Check for an interrupt, if applicable
    if (mIsInterruptable) {
#ifdef CTRLC
      // CTRL-C method
      if (utIsInterruptPending()) {
	mTaskMutex.lock();  // lock so mIsRunning doesn't change in middle of task-retrieval
	mIsRunning = false;
	mTaskMutex.unlock();
      }
#else
      // Dialog method: see if the dialog is still there
      // Because Matlab is not thread-safe, run this only from the main thread
      if (boost::this_thread::get_id() == mMainThreadID) {
	const mxArray *ret = mexGet(mxGetScalar(mInterruptHandle),"Parent");
	if (ret == NULL) {
	  mTaskMutex.lock();
	  mIsRunning = false;
	  mTaskMutex.unlock();
	}
      }
#endif  // CTRLC
    }
#endif  // mex_h
    // Acquire the task lock
    mTaskMutex.lock();
    // Check to see if the queue is empty
    if (mTaskQueue.empty())
      mIsRunning = false;
    bool retval = mIsRunning;  // make a thread-local copy of the running status
    int nTasks = 0;
    if (mIsRunning) {
      // Retrieve thread-local copy of the # of tasks remaining
      // (used for showing progress)
      nTasks = mTaskQueue.size();
      // Deliver the next task
      T = mTaskQueue.front();
      mTaskQueue.pop();
    }
    // Release the task lock
    mTaskMutex.unlock();
    // Produce the console output
    // To be thread-safe, we run this only in the main thread
    if (boost::this_thread::get_id() == mMainThreadID) {
      // Display any accumulated messages
      mMessageMutex.lock();    // acquire the message lock
      while (!mMessageQueue.empty()) {
#ifdef mex_h
	mexPrintf("%s\n",mMessageQueue.front().c_str());
#else
	std::cout << mMessageQueue.front() << std::endl;
#endif
	mMessageQueue.pop();
      }
      mMessageMutex.unlock();
      // Show progress, if applicable
      if (mIsShowProgress && mLastReport-nTasks >= mProgressInterval) {
	mLastReport = nTasks;
#ifdef mex_h
	mexPrintf("%d...",nTasks);
	mexEvalString("drawnow;");  // flush screen output
#else
	std::cout << nTasks << "..." << std::flush;
#endif
      }
    }
    return retval;  // return the running status
  }
  void queueMessage(const std::string &s) {
    mMessageMutex.lock();
    mMessageQueue.push(s);
    mMessageMutex.unlock();
  }
  int progressInterval() const { return mProgressInterval; }
  void progressInterval(int i) { mProgressInterval = i; }
  int size() const { return mTaskQueue.size(); }
  bool isMainThread() const { return boost::this_thread::get_id() == mMainThreadID; }

private:
  std::queue<TaskType>    mTaskQueue;
  boost::mutex            mTaskMutex;
  std::queue<std::string> mMessageQueue;
  boost::mutex            mMessageMutex;
  bool                    mIsRunning;
  bool                    mIsShowProgress;
  bool                    mIsInterruptable;
  int                     mProgressInterval;
  int                     mLastReport;
  boost::thread::id       mMainThreadID;
#if defined(mex_h) && !defined(CTRLC)
  mxArray*                mInterruptHandle;
#endif
};

#endif // BOOSTTHREADPOOL_H
