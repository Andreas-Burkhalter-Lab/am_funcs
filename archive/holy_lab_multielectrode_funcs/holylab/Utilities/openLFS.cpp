/*=================================================================
  open file with LFS mode

  todo: currently only read file
 *=================================================================*/

#include "../commong/lfs_common.h"
#include "errno.h"

#if defined(WIN32) || defined(__WIN32) || defined(WIN)
#   define _sys_errlist sys_errlist
#endif

#include "mex.h"

#include "matlab_arg.h"

//---------------------------------------------------------------------

int openLFS(const char* filename)
{
   int fd=open(filename, O_RDONLY| O_LARGEFILE);

   return fd;
}

//---------------------------------------------------------------------
void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* Check for proper number of input and  output arguments */    
    if (nrhs !=1) {
        mexErrMsgTxt("one input argument required.");
    } 
    if(nlhs > 2){
        mexErrMsgTxt("0~2 output arguments required.");
    }

    /* input must be a string */
    if( mxIsChar(prhs[0]) != 1)  mexErrMsgTxt("Input must be a string.");

    /* input must be a row vector */
    if(mxGetM(prhs[0])!=1)  mexErrMsgTxt("Input must be a row vector.");
    
    /* get the length of the input string */
    int buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

    /* allocate memory for input and output strings */
    char * input_buf=(char*)mxCalloc(buflen, sizeof(char));
    
    /* copy the string data from prhs[0] into a C string input_ buf.
     * If the string array contains several rows, they are copied,
     * one column at a time, into one long string array.
     */
    int status = mxGetString(prhs[0], input_buf, buflen);
    if(status != 0) mexWarnMsgTxt("Not enough space. String is truncated.");
    
    int fd=openLFS(input_buf);
    if(fd<0){mexErrMsgTxt("error when open the file");}
  
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = fd;

    if(nlhs==2){
       string errmsg=sys_errlist[errno];
       matlab_arg::setStringArg(plhs[1], errmsg);
    }

}

