#include "mex.h"

#include <iostream>

using namespace std;


#include <sys/types.h>
#include <unistd.h>

/*
fid: file identifier, returned by matlab fopen()
fd: file descriptor, returned by c/c++ open()

purpose:
test if fid returned by matlab fopen() can be treated as fd and be used by read() in c/c++ code

Answer: no, unfortunately
 */

void readfirstline(int fd)
{
   char buf[1024];
   off_t posNow=lseek(fd, 0, SEEK_SET);
   int nRead=read(fd, buf, sizeof(buf));
   if(nRead<1){mexErrMsgTxt("read error"); }
   buf[nRead-1]=0;
   mexPrintf("buf is:\n%s\n",buf);
   
}


//input arg: a fd can be used by c code

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x,*y;
  int     mrows,ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=1) {
    mexErrMsgTxt("One input required.");
  } 
  if(nlhs>0) {
    mexErrMsgTxt("no output arguments is required");
  }
  
  /* The input must be a noncomplex scalar double.*/
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
  
  /* Assign pointers to each input and output. */
  int fd = *mxGetPr(prhs[0]);
  mexPrintf("fd=%d\n",fd);
  
  /* then read the first line in file  */
  readfirstline(fd);
}
