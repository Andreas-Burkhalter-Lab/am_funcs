#include "mex.h"
#include <string.h>

const int mode_full = 0;
const int mode_diag = 1;

/*
 * This is the Matlab wrapper
 */

void gradwork(double *x,double *y,double *n,double *G,double *Qcv,int N,int d,int dp,double *grad,int mode);

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *x, *y, *G, *Qcv, *n, *grad;
  int N,d,dp,mode;
  const mxArray *curarg;
  const int bufsize = 10;
  char modestr[bufsize];

  if (nrhs < 5)
    mexErrMsgTxt("gradwllcv: requires five inputs");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("gradwllcv: x must be a real matrix");
  x = mxGetPr(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("gradwllcv: y must be a real matrix");
  y = mxGetPr(curarg);
  dp = mxGetM(curarg);
  if (N != mxGetN(curarg))
    mexErrMsgTxt("gradwllcv: the number of points in x and y disagree");

  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("gradwllcv: n must be a real vector");
  n = mxGetPr(curarg);
  if (mxGetM(curarg)*mxGetN(curarg) != N) {
    mexPrintf("%d %d\n",mxGetM(curarg),mxGetN(curarg));
    mexErrMsgTxt("gradwllcv: n must have the same number of points as x & y");
  }

  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("gradwllcv: G must be a real matrix");
  G = mxGetPr(curarg);
  if (mxGetN(curarg) != N || mxGetN(curarg) != N)
    mexErrMsgTxt("gradwllcv: G must be N-by-N");

  curarg = prhs[4];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("gradwllcv: Qcv must be a real vector");
  Qcv = mxGetPr(curarg);
  if (mxGetM(curarg)*mxGetN(curarg) != N)
    mexErrMsgTxt("gradwllcv: Qcv must have the same number of points as x & y");

  curarg = prhs[5];
  if (!mxIsChar(curarg))
    mexErrMsgTxt("gradwllcv: mode must be a string");
  mxGetString(curarg,modestr,bufsize);
  if (strncmp(modestr,"full",bufsize) == 0)
    mode = mode_full;
  else if (strncmp(modestr,"diag",bufsize) == 0)
    mode = mode_diag;
  else
    mexErrMsgTxt("gradwllcv: mode not recognized");

  //mexPrintf("d %d, N %d, q %d, nlhs %d, nrhs %d\n",d,N,q,nlhs,nrhs);

  // Set up the output
  if (mode == mode_full)
    plhs[0] = mxCreateDoubleMatrix(dp,d,mxREAL);
  else if (mode == mode_diag)
    plhs[0] = mxCreateDoubleMatrix(dp,1,mxREAL);
  grad = mxGetPr(plhs[0]);
  // Do the actual work
  gradwork(x,y,n,G,Qcv,N,d,dp,grad,mode);
  return;
}


/*
 * 
 * d-by-N, y is d-by-q, and both are organized by column in memory.
 */
void gradwork(double *x,double *y,double *n,double *G,double *Qcv,int Nu,int d,int dp,double *grad,int mode)
{
  double gradtmp;
  double *xp1,*yp1,*xp2,*yp2,*Gp;
  int i,j,N;
  int kap1,kap2;

  // Calculate total number of points, including multiplicities
  N = 0;
  for (i = 0; i < Nu; i++)
    N = N+int(n[i]);

  if (mode == mode_full) {
    // Loop over the coordinates
    for (kap1 = 0; kap1 < dp; kap1++) {  // First index of kappa
      for (kap2 = 0; kap2 < d; kap2++) { // Second index of kappa
	gradtmp = 0;
	xp1 = x + kap2;  // offset to pick correct coordinate
	yp1 = y + kap1;  // offset to pick correct coordinate
	Gp = G;
	// Loop over the point pairs
	for (i = 0; i < Nu; i++,xp1+=d,yp1+=dp) {
	  //xp2 = xp1;
	  //yp2 = yp1;
	  xp2 = x + kap2;
	  yp2 = y + kap1;
	  for (j = 0; j < Nu; j++,xp2+=d,yp2+=dp,Gp++) {
	    gradtmp += (*xp2-*xp1) * *Gp * (*yp2-*yp1) * n[i] * n[j] / Qcv[i];
	  }
	}
	grad[kap1+kap2*dp] = gradtmp/N;
      }
    }
  }
  else if (mode == mode_diag) {
    // Loop over the coordinates
    for (kap1 = 0; kap1 < dp; kap1++) {  // First index of kappa
      gradtmp = 0;
      xp1 = x + kap1;  // offset to pick correct coordinate
      yp1 = y + kap1;  // offset to pick correct coordinate
      Gp = G;
      // Loop over the point pairs
      for (i = 0; i < Nu; i++,xp1+=d,yp1+=dp) {
	//xp2 = xp1;
	//yp2 = yp1;
	xp2 = x + kap1;
	yp2 = y + kap1;
	for (j = 0; j < Nu; j++,xp2+=d,yp2+=dp,Gp++) {
	  gradtmp += (*xp2-*xp1) * *Gp * (*yp2-*yp1) * n[i] * n[j] / Qcv[i];
	}
      }
      grad[kap1] = gradtmp/N;
    }
  }
}
