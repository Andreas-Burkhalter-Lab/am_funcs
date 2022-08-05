#include "mex.h"

/*
 * This is the Matlab wrapper
 */

void meanshift1dloop_work(double *x,double *xo,int N,double thresh);

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *x, *xo, thresh;
  int N;
  const mxArray *curarg;

  if (nrhs != 2)
    mexErrMsgTxt("meanshift1d1: requires two inputs");
  if (nlhs != 1)
    mexErrMsgTxt("meanshift1d1: requires one output");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) ||
      mxGetNumberOfDimensions(curarg) > 2 ||
      (mxGetM(curarg) != 1 && mxGetN(curarg) != 1))
    mexErrMsgTxt("meanshift1d1: x must be a real vector");
  x = mxGetPr(curarg);
  N = mxGetNumberOfElements(curarg);

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) ||
      mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("meanshift1d1: thresh must be a real scalar");
  thresh = mxGetScalar(curarg);

  // Set up the output
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL);
  xo = mxGetPr(plhs[0]);

  // Do the actual work
  meanshift1dloop_work(x,xo,N,thresh);
  return;
}


/*
 * Given a (sorted) list of points and a threshold, output a new set
 * of points, where each point is replaced by the center of mass of
 * all points within a given distance (thresh).
 *
 */
void meanshift1dloop_work(double *x,double *xo,int N,double thresh)
{
  double *xc, *xl, *xr, *xend, *xendm1;
  xend = x + N;
  xendm1 = xend-1;

  for (xc = x; xc < xend; xc++,xo++) {
    xl = xr = xc;
    *xo = *xc;
    // Move right boundary as far as possible
    while (xr < xendm1 && *(xr+1) - *xc < thresh) {
      xr++;
      *xo += *xr;
    }
    // Move left boundary as far as possible
    while (xl > x && *xc - *(xl-1) < thresh) {
      xl--;
      *xo += *xl;
    }
    *xo = *xo/(xr-xl+1);
  }
}
