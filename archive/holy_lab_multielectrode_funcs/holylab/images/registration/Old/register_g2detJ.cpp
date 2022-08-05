#include "mex.h"
#include <math.h>
#include "image_utilities.h"

/*
 * register_g2detJ: calculate detJ without storing J
 *
 * Syntax:
 *   detJ = register_g2detJ(g1,g2,...)
 *   detJ = register_g2detJ(g1,g2,...,handle_edges)
 * where
 *   g1, g2 are the deformation arrays, i.e., g{:}.
 *   handle_edges is a boolean (default true) that, if true, results
 *     in finite values being produced on edges (if false, NaNs);
 * and
 *   detJ is the determinant of the Jacobian at each spatial location.
 *
 * Copyright 2006-2007 by Timothy E. Holy
 */

void g2detJwork(int n_dims,const float *g[],const int *sz_g,float *detJ,bool handle_edges);

/*
 * This is the Matlab wrapper
 */

#define MAX_DIMS 3

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *g[MAX_DIMS];
  float *detJ;
  int n_dims,n_dims_tmp,dimIndex;
  int sz_g[MAX_DIMS],sz_g_tmp[MAX_DIMS];
  int *sz_g_ptr;
  const mxArray *curarg;
  bool handle_edges;

  if (nrhs < 1)
    mexErrMsgTxt("register_g2detJ: requires at least one input");
  if (nlhs != 1)
    mexErrMsgTxt("register_g2detJ: requires one output");

  // Parse the inputs
  // handle_edges: check whether the last item is a scalar
  handle_edges = true;
  curarg = prhs[nrhs-1];
  if (mxGetNumberOfElements(curarg) == 1) {
    handle_edges = (mxGetScalar(curarg) > 0);
    n_dims = nrhs-1;
  }
  else
    n_dims = nrhs;

  // Components of g
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    curarg = prhs[dimIndex];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
      mexErrMsgTxt("register_g2detJ: components of g must be a real single-precision array");
    g[dimIndex] = (float *) mxGetData(curarg);
    if (dimIndex == 0)
      sz_g_ptr = sz_g;
    else
      sz_g_ptr = sz_g_tmp;
    n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				       mxGetNumberOfDimensions(curarg),
				       sz_g_ptr,
				       MAX_DIMS);
    if (n_dims_tmp != n_dims)
      mexErrMsgTxt("register_g2detJ: dimensionality of components of g must equal image dimensionality");
    if (dimIndex != 0)
      if (!validate_dimensions(sz_g,n_dims,sz_g_tmp,n_dims_tmp))
	mexErrMsgTxt("register_g2detJ: all components of g must have the same size");
  }
    
  // Set up the output. Make it the same geometry as the input,
  // including singleton dimensions.
  plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
				 mxGetDimensions(prhs[0]),
				 mxSINGLE_CLASS,mxREAL);
  detJ = (float *) mxGetData(plhs[0]);

  // Do the actual work
  g2detJwork(n_dims,g,sz_g,detJ,handle_edges);

  return;
}


/*
 * Do the actual interpolation.
 */
void g2detJwork(int n_dims,const float *g[],const int *sz_g,float *detJ,bool handle_edges)
{
  long n_pixels_out,pixel_step[3];
  int dimIndex1,dimIndex2,coord[3];
  float *detJ_end, J[3][3];
  bool is_nan;

  n_pixels_out = 1;
  for (dimIndex1 = 0; dimIndex1 < n_dims; dimIndex1++)
    n_pixels_out *= sz_g[dimIndex1];
  detJ_end = detJ+n_pixels_out;

  pixel_step[0] = 1;
  for (dimIndex1 = 1; dimIndex1 < n_dims; dimIndex1++)
    pixel_step[dimIndex1] = pixel_step[dimIndex1-1] * sz_g[dimIndex1-1];

  // Initialize coords
  for (dimIndex1 = 0; dimIndex1 < n_dims; dimIndex1++)
    coord[dimIndex1] = 0;

  for (; detJ < detJ_end; detJ++) {
    // Calculate the Jacobian at this point. Use centered differencing
    // where possible, one-sided differencing at the edges.
    is_nan = false;
    for (dimIndex1 = 0; dimIndex1 < n_dims; dimIndex1++) {
      for (dimIndex2 = 0; dimIndex2 < n_dims; dimIndex2++) {
	if (coord[dimIndex2] > 0 && coord[dimIndex2] < sz_g[dimIndex2]-1) {
	  J[dimIndex1][dimIndex2] = (*(g[dimIndex1] + pixel_step[dimIndex2]) -
				     *(g[dimIndex1] - pixel_step[dimIndex2]))/2;
	}
	else if (!handle_edges)
	  is_nan = true;
	else if (coord[dimIndex2] == 0)
	  J[dimIndex1][dimIndex2] = *(g[dimIndex1] + pixel_step[dimIndex2]) -
	    *(g[dimIndex1]);
	else
	  J[dimIndex1][dimIndex2] = *(g[dimIndex1]) -
	    *(g[dimIndex1] - pixel_step[dimIndex2]);
      }
    }
    // Calculate the determinant
    if (is_nan)
      *detJ = mxGetNaN();
    else if (n_dims == 1)
      *detJ = J[0][0];
    else if (n_dims == 2)
      *detJ = J[0][0] * J[1][1] - J[1][0] * J[0][1];
    else
      *detJ = J[0][0]*(J[1][1]*J[2][2] - J[2][1]*J[1][2]) -
	J[0][1]*(J[1][0]*J[2][2] - J[2][0]*J[1][2]) +
	J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);
    // Update the deformation pointers
    for (dimIndex1 = 0; dimIndex1 < n_dims; dimIndex1++)
      (g[dimIndex1])++;
    // Update the coordinate counters
    dimIndex1 = 0;
    while (dimIndex1 < n_dims) {
      if (++coord[dimIndex1] >= sz_g[dimIndex1]) {
	coord[dimIndex1] = 0;
	dimIndex1++;
      }
      else
	break;
    }
  }
}

