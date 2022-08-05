/*=================================================================

  close file opened with LFS mode

 *=================================================================*/

#include "../commong/lfs_common.h"

#include "mex.h"

//---------------------------------------------------------------------

int closeLFS(int fd)
{
   return close(fd);
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
  
   // get the file descriptor:
   int fd = (int) *mxGetPr(prhs[0]); //todo: may lose precision?
  
  /* Call the timestwo subroutine. */
   closeLFS(fd);

}

