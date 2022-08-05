# include "mex.h"
# include <iostream>
# include "matrix.h"
# include <sys/time.h>
# include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sched.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <signal.h>
#include <errno.h>

#define DEBUG 1

static double *valveTimeSeq = 0;
static int ncols;
static int tCounter;
static int vCounter;

/*
 * This function gets called whenever the timer is up
 */
void my_alrm_handler(int num)
{
  struct itimerval itimer,oldtimer;
  
  if (DEBUG) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    mexPrintf("Alarm response time: %d\n",tv.tv_sec+double(tv.tv_usec)/1e6);
  }
  
  long sec, microsec;
  sec = 0;
  microsec = 0;
  
  tCounter = tCounter + 2;
  
  //mexPrintf("ncols is %d\n",ncols);
  
  if (pr[vCounter]== 0){
    mexPrintf("Open valve %f\n",valveTimeSeq[vCounter+2]);
  }
  else {
    mexPrintf("Close valve %f\n",valveTimeSeq[vCounter]);
	}

	mexPrintf("vCounter is %d\n",vCounter);
	vCounter = vCounter + 2;

	if (tCounter < 2*ncols){
	   mexPrintf("Wait is %f\n",valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2]);
	   sec = long (valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2]);
	   microsec = long(((valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2])-sec)*1000000);
	   mexPrintf("sec is %d\n",sec);
	   mexPrintf("microsec is %d\n",microsec);
	   itimer.it_interval.tv_sec = 0;
  	   itimer.it_interval.tv_usec= 0;
	   itimer.it_value.tv_sec = sec;
  	   itimer.it_value.tv_usec = microsec;
	   setitimer(ITIMER_REAL, &itimer, &oldtimer);
	   }
	return;
	// If we detect we're at the end of the sequence, free the memory at valveTimeSeq.  Set it to 0 to indicate it's free.
}


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

  int nelements;
  double *pr;
  int mrows, i;
  
  tCounter = 1;
  vCounter = 0;
  
  // If there is a job already running, we must first terminate it
  if (valveTimeSeq) {
    do something here
  }

  // If there are no arguments, the user just wanted to terminate.
  if (nrhs < 1)
    return;

  // The input must be a 2xn double real matrix
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (!mxIsDouble(prhs[0]) || (mrows != 2) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("Input must be a 2xn double real matrix");
  
  // Get the data, and store it in a new matrix, to make sure this
  // memory doesn't get de-allocated in MATLAB
  pr = (double *)mxGetPr(prhs[0]);
  nelements = mrows*ncols;
  valveTimeSeq = (double *) malloc(nelements*sizeof(double));
  for (i = 0; i < nelements; i++)
    valveTimeSeq[i] = pr[i];
  
  
  // Set up the job scheduling
  struct itimerval itimer,oldtimer;
  struct sigaction alrm_action,old_action;
  
  // Initialize the sigaction structure; make sure that
  // the only signal it responds to is when the alarm timer goes off
  sigemptyset(&alrm_action.sa_mask);
  sigaddset(&(alrm_action.sa_mask), SIGALRM);
  // Automatically restart any system calls interrupted by
  // this signal
  alrm_action.sa_flags = SA_RESTART;
  // Specify the function that should be called when the alarm goes off
  alrm_action.sa_handler = my_alrm_handler;
  // Parameters are set. Add this to the list of signals this process
  // can handle.
  if (sigaction(SIGALRM, &alrm_action, &old_action) < 0)
    mexErrMsgTxt("Error registering the signal");
  // Make sure we aren't blocking this signal
  sigprocmask(SIG_UNBLOCK, &alrm_action.sa_mask, &old_action.sa_mask);

	//#define printvar(x)  printf(#x "= %f\n", (double) x)

	tCounter = tCounter + 2;

	long sec =0;
	long microsec = 0;

	sec = long (valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2]);
	microsec = long(((valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2])-sec)*1000000);


	mexPrintf("sec is %d\n",sec);
	mexPrintf("microsec is %d\n",microsec);

	mexPrintf("Wait is %f\n",valveTimeSeq[tCounter]-valveTimeSeq[tCounter-2]);
  	itimer.it_value.tv_sec = sec;
  	itimer.it_value.tv_usec = microsec;
	itimer.it_interval.tv_sec = 0;
  	itimer.it_interval.tv_usec= 0;

	if (DEBUG){
	   struct timeval tv;
	   struct timezone tz;
	   gettimeofday(&tv,&tz);
	   printvar(tv.tv_sec);
 	   printvar(tv.tv_usec);
	   }

	setitimer(ITIMER_REAL, &itimer, &oldtimer);

	printvar(itimer.it_interval.tv_sec);
  	printvar(itimer.it_interval.tv_usec);
	printvar(itimer.it_value.tv_sec);
  	printvar(itimer.it_value.tv_usec);
}
