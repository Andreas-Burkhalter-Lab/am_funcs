#include "mex.h"

/*
 * register_detJ: calculate determinant arrays
 *
 * Syntax:
 *   detJ = register_detJ(J)
 * where
 *   J is a n-by-n-by-[sz] array
 * calculates
 *   detJ, an array of size [sz], where each element of detJ, say detJ(i),
 *     is the determinant of the matrix J(:,:,i).
 *
 * Copyright 2006 by Timothy E. Holy
 */

void detJwork(float *J,float *Jend,int matrix_dim,float *detJ);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  float *J,*detJ;
  int n_dims,matrix_dim,n_spatial_dims,spatial_sz[3],i;
  const int *sz;
  const mxArray *curarg;

  if (nrhs != 1)
    mexErrMsgTxt("register_detJ: requires one input");
  if (nlhs != 1)
    mexErrMsgTxt("register_detJ: requires one output");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_detJ: J must be a real single-precision array");
  J = (float *) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  if (n_dims < 2)
    mexErrMsgTxt("register_detJ: J must have at least two dimensions");
  sz = mxGetDimensions(curarg);

  // Sanity check on dimensionality
  matrix_dim = sz[0];
  if (matrix_dim != sz[1])
    mexErrMsgTxt("register_detJ: J must be a n-by-n-by-[sz] array");
  if (matrix_dim > 3)
    mexErrMsgTxt("register_detJ: for now, handles up to 3-by-3 matrices only");

  // Set up the outputs
  n_spatial_dims = n_dims-2;
  for (i = 0; i < n_spatial_dims; i++)
    spatial_sz[i] = sz[2+i];
  if (n_spatial_dims < 2) {
    for (i = n_spatial_dims; i < 3; i++)
      spatial_sz[i] = 1;
    n_spatial_dims = 2;
  }

  plhs[0] = mxCreateNumericArray(n_spatial_dims,spatial_sz,mxSINGLE_CLASS,mxREAL);
  detJ = (float *) mxGetData(plhs[0]);

  // Do the actual work
  detJwork(J,J+mxGetNumberOfElements(prhs[0]),matrix_dim,detJ);

  return;
}


/*
 * Do the actual computation.
 * It uses explicit formulas for determinants in dimensions 1, 2, and 3. 
 * This is a big speed boost.
 */
void detJwork(float *J,float *Jend,int matrix_dim,float *detJ)
{
  int pixelSkip;

  pixelSkip = matrix_dim*matrix_dim;
  if (matrix_dim == 1)
    for (; J < Jend; J++,detJ++)
      *detJ = *J;
  else if (matrix_dim == 2)
    for (; J < Jend; J += pixelSkip, detJ++)
      *detJ = J[0]*J[3]-J[1]*J[2];
  else if (matrix_dim == 3)
    for (; J < Jend; J += pixelSkip, detJ++)
      *detJ = J[0]*(J[4]*J[8] - J[5]*J[7]) - J[3]*(J[1]*J[8] - J[2]*J[7])
	+ J[6]*(J[1]*J[5] - J[2]*J[4]);
  return;
}
