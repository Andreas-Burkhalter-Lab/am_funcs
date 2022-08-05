//compile:
//linux:      mex sleep_g.cpp
//window:     mex -DWIN sleep_g.cpp

#include "mex.h"

#if defined(WIN32) || defined(__WIN32) || defined(WIN)
#   include <windows.h>
#else
#   include <time.h>
#endif
//---------------------------------------------------------------------

static void sleep_g(double howlong)
{
#if defined(WIN32) || defined(__WIN32) || defined(WIN)
   Sleep(howlong*1000);	
#else //else: linux
   struct timespec t;
   t.tv_sec=howlong;
   t.tv_nsec=(howlong-t.tv_sec)*1e9;
   while(nanosleep(&t, &t) ){
      //nop
   }
#endif   
}

//---------------------------------------------------------------------
void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if(nrhs!=1) {
      mexErrMsgTxt("One input required.");
   } 
   if(nlhs>0) {
      mexErrMsgTxt("no output arguments is required");
   }
  
   /* The input must be a noncomplex scalar double.*/
   int mrows = mxGetM(prhs[0]);
   int ncols = mxGetN(prhs[0]);
   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
       !(mrows==1 && ncols==1) ) {
      mexErrMsgTxt("Input must be a noncomplex scalar double.");
   }
  
   double duration = *mxGetPr(prhs[0]); 
  
   sleep_g(duration);

}

