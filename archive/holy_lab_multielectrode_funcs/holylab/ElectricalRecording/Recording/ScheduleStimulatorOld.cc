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

static double *pr;
static int ncols;
static int tCounter;
static int vCounter;

/*
 * This function gets called whenever the timer is up
 */
void my_alrm_handler(int num)
{
#define printvar(x)  printf(#x "= %f\n", (double) x)
	struct itimerval itimer,oldtimer;

	if (DEBUG) {
	  struct timeval tv;
	  struct timezone tz;
	  gettimeofday(&tv,&tz);
	  printvar(tv.tv_sec);
	  printvar(tv.tv_usec);
	  mexPrintf("We got to the alarm handler\n");
	}

	long sec, microsec;
	sec = 0;
	microsec = 0;

	tCounter = tCounter + 2;

	//mexPrintf("ncols is %d\n",ncols);

	if (pr[vCounter]== 0){
	   mexPrintf("Open valve %f\n",pr[vCounter+2]);
	}
	else {
	   mexPrintf("Close valve %f\n",pr[vCounter]);
	}

	mexPrintf("vCounter is %d\n",vCounter);
	vCounter = vCounter + 2;

	if (tCounter < 2*ncols){
	   mexPrintf("Wait is %f\n",pr[tCounter]-pr[tCounter-2]);
	   sec = long (pr[tCounter]-pr[tCounter-2]);
	   microsec = long(((pr[tCounter]-pr[tCounter-2])-sec)*1000000);
	   mexPrintf("sec is %d\n",sec);
	   mexPrintf("microsec is %d\n",microsec);
	   itimer.it_interval.tv_sec = 0;
  	   itimer.it_interval.tv_usec= 0;
	   itimer.it_value.tv_sec = sec;
  	   itimer.it_value.tv_usec = microsec;
	   setitimer(ITIMER_REAL, &itimer, &oldtimer);
	   }
	return;
}


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

	int elements;
	double *pi,*pind;
	int mrows, i;

	tCounter = 1;
	vCounter = 0;

	// Check for proper number of arguments
	if (nrhs != 1)
	   mexErrMsgTxt("One input required");
	if (nlhs > 1)
	   mexErrMsgTxt("Too many output arguments");

	// The input must be a 2xn double matrix
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (!mxIsDouble(prhs[0])|| (mrows != 2))
	   mexErrMsgTxt("Input must be a 2xn double matrix");

	mexPrintf("We got to spot 1\n");

	// Get the data
	pr = (double *)mxGetPr(prhs[0]);
	pi = (double *)mxGetPi(prhs[0]);

	// Allocate space for the return argument
	plhs[0] = mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	mexPrintf("We got to spot 2\n");
	pind = mxGetPr(plhs[0]);

	elements = mxGetNumberOfElements(prhs[0]);

	for(i=0;i<elements;i++){
	   pind[i] = pr[i];
	   }
	mexPrintf("We got to spot 3\n");

	mexPrintf("We got to spot 4\n");


	struct itimerval itimer,oldtimer;
	struct sigaction alrm_action,old_action;

	sigemptyset(&alrm_action.sa_mask);
	sigaddset(&(alrm_action.sa_mask), SIGALRM);
	alrm_action.sa_flags = SA_RESTART;
	alrm_action.sa_handler = my_alrm_handler;
	sigaction(SIGALRM, &alrm_action, &old_action);
	sigprocmask(SIG_UNBLOCK, &alrm_action.sa_mask, &old_action.sa_mask);

#define printvar(x)  printf(#x "= %f\n", (double) x)

	tCounter = tCounter + 2;

	long sec =0;
	long microsec = 0;

	sec = long (pr[tCounter]-pr[tCounter-2]);
	microsec = long(((pr[tCounter]-pr[tCounter-2])-sec)*1000000);


	mexPrintf("sec is %d\n",sec);
	mexPrintf("microsec is %d\n",microsec);

	mexPrintf("Wait is %f\n",pr[tCounter]-pr[tCounter-2]);
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
