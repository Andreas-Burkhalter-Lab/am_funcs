#include "mex.h"

/* adaptive_meanshift_step: take one step of an adaptive meanshift
 *
 * Syntax:
 *   [yout,R2,n] = adaptive_meanshift_step(x,y,factor,minN)
 * where
 *   x is a d-by-N matrix of data points in d-dimensional space;
 *   y is a d-by-q matrix of landmarks that will move up the density
 *      defined by x;
 *   factor (default 3) is the number of times that the step size needs
 *      to exceed the standard error by to be considered significant;
 *   minN (default 10) is the minimum number of neighbors to use to
 *      initially judge whether a step should be taken.
 * and
 *   yout is the new set of landmark positions;
 *   R2 is a 1-by-q vector containing the set of square radii used
 *     for each landmark;
 *   n is a 1-by-q vector containing the number of data points (x)
 *     within the radius at each landmark.

/*
 * This is the Matlab wrapper
 */

int is_scalar(const mxArray *m);
void heapsort_index(double *val,long *index,long n);
void amswork(const double *x,const double *y,int N,int d,int q,double factor2,int minN,double *yout,double *R2,double *n,double *deltax,double *sqrdist,long *nbrIndex,double *deltax_accum,double *stepsq_accum,double *sumsq_accum);

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *x, *y, factor, *yout, *R2, *n;
  double *deltax, *sqrdist, *deltax_accum, *sumsq_accum, *stepsq_accum;
  long *nbrIndex;
  int N,d,q,minN;
  const mxArray *curarg;

  if (nrhs < 2)
    mexErrMsgTxt("adaptive_meanshift_step: requires at least two inputs");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("adaptive_meanshift_step: x must be a real matrix");
  x = mxGetPr(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("adaptive_meanshift_step: y must be a real matrix");
  y = mxGetPr(curarg);
  if (d != mxGetM(curarg))
    mexErrMsgTxt("adaptive_meanshift_step: the dimensionality of x and y disagree");
  q = mxGetN(curarg);

  if (nrhs > 2) {
    curarg = prhs[2];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !is_scalar(curarg))
      mexErrMsgTxt("adaptive_meanshift_step: factor must be a real scalar");
    factor = mxGetScalar(curarg);
  } else
    factor = 3.0;
    
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !is_scalar(curarg))
      mexErrMsgTxt("adaptive_meanshift_step: minN must be a real scalar");
    minN = int(mxGetScalar(curarg));
    if (minN < 2)
      mexErrMsgTxt("adaptive_meanshift_step: minN must be >= 2");
    if (minN > N)
      mexErrMsgTxt("adaptive_meanshift_step: minN is larger than N");
  } else
    minN = (N >= 10) ? 10 : N;  // set to min(10,N)

  //mexPrintf("d %d, N %d, q %d, nlhs %d, nrhs %d\n",d,N,q,nlhs,nrhs);

  // Set up the outputs
  plhs[0] = mxCreateDoubleMatrix(d,q,mxREAL);
  yout = mxGetPr(plhs[0]);
  R2 = NULL;
  n = NULL;
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,q,mxREAL);
    R2 = mxGetPr(plhs[1]);
  } else
    R2 = NULL;
  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(1,q,mxREAL);
    n = mxGetPr(plhs[2]);
  } else
    n = NULL;

  // Allocate temporary storage
  deltax = (double *) mxMalloc(d*N*sizeof(double));
  sqrdist = (double *) mxMalloc(N*sizeof(double));
  nbrIndex = (long *) mxMalloc(N*sizeof(long));
  deltax_accum = (double *) mxMalloc(d*N*sizeof(double));
  sumsq_accum = (double *) mxMalloc(N*sizeof(double));
  stepsq_accum = (double *) mxMalloc(N*sizeof(double));

  // Do the actual work
  amswork(x,y,N,d,q,factor*factor,minN,yout,R2,n,deltax,sqrdist,nbrIndex,deltax_accum,stepsq_accum,sumsq_accum);

  // Free temporary storage
  mxFree(deltax);
  mxFree(sqrdist);
  mxFree(nbrIndex);
  mxFree(deltax_accum);
  mxFree(sumsq_accum);
  mxFree(stepsq_accum);

  return;
}

int is_scalar(const mxArray *m)
{
  return (mxGetNumberOfElements(m) == 1);
}

void amswork(const double *x,const double *y,int N,int d,int q,double factor2,int minN,double *yout,double *R2,double *n,double *deltax,double *sqrdist,long *nbrIndex,double *deltax_accum,double *stepsq_accum,double *sumsq_accum)
{
  const double *xend,*yend,*xp;
  double *deltaxp,*sqrdistp,*deltax_accump;
  long xIndex,dimIndex;

  yend = y + d*q;
  xend = x + d*N;
  for ( ; y < yend; y += d, yout += d) {
    // y points to the current landmark
    // Start x at the beginning of the data points on each iteration
    // Calculate deltax, and the square distance between the current landmark and the current point
    for (xp = x, deltaxp = deltax, sqrdistp = sqrdist; xp < xend; sqrdistp++) {
      *sqrdistp = 0;
      for (dimIndex = 0; dimIndex < d; dimIndex++, xp++, deltaxp++) {
	*deltaxp = *xp - y[dimIndex];
	*sqrdistp += *deltaxp * *deltaxp;
      }
      //mexPrintf("  %g\n",*sqrdistp);
    }
    // Sort the data points in order of increasing square distance
    heapsort_index(sqrdist,nbrIndex,N);
    // Now march through the data, comparing whether
    //   (sum(displacement))^2   >   sum(displacement^2)
    // by a sufficient factor.
    for (dimIndex = 0; dimIndex < d; dimIndex++) {
      // initialize with the closest point
      deltax_accum[dimIndex] = deltax[d*nbrIndex[0]+dimIndex];
      sumsq_accum[0] = sqrdist[0];
      stepsq_accum[0] = sqrdist[0];
    }
    for (xIndex = 1, deltax_accump = deltax_accum+d; xIndex < N; xIndex++) {
      stepsq_accum[xIndex] = 0;
      for (dimIndex = 0; dimIndex < d; dimIndex++,deltax_accump++) {
	*deltax_accump = *(deltax_accump - d) + 
	  deltax[d*nbrIndex[xIndex] + dimIndex];
	stepsq_accum[xIndex] += *deltax_accump * *deltax_accump;
      }
      sumsq_accum[xIndex] = sumsq_accum[xIndex-1] + 
	sqrdist[xIndex];
      //mexPrintf("  nbr %d: sqrdist %g, sumsq %g, stepsq %g\n",nbrIndex[xIndex],sqrdist[xIndex],sumsq_accum[xIndex],stepsq_accum[xIndex]);
      // Are we above threshold?
      if (xIndex >= minN) // do this only after minN points
	if (stepsq_accum[xIndex] > factor2*sumsq_accum[xIndex])
	  break;
    }
    // If we hit threshold, then march backwards until the comparison
    // just creeps above threshold (with unity factor). Otherwise, stay
    // at the last point.
    //mexPrintf("xIndex: %d\n",xIndex);
    if (xIndex < N)
      while (xIndex > 0 && stepsq_accum[xIndex-1] > sumsq_accum[xIndex-1])
	xIndex--;
    else
      xIndex = N-1;
    //mexPrintf("xIndex now: %d\n",xIndex);
    //mexPrintf("dxaccum: %g, y: %g\n",deltax_accum[xIndex],y[0]);
    // OK. Report the answer: first, the new y position. We have to
    // add the current position to the displacement
    for (dimIndex = 0; dimIndex < d; dimIndex++)
      yout[dimIndex] = deltax_accum[d*xIndex+dimIndex]/(xIndex+1) + y[dimIndex];
    // Now report the # of points and the radius squared
    if (n != NULL) {
      *n = xIndex+1;
      n++;
    }
    if (R2 != NULL) {
      *R2 = sqrdist[xIndex];
      R2++;
    }
  } // End loop over landmarks
}
