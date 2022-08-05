#include "mex.h"

/*
 * register_Jtinvw: solve arrays of linear equations
 *
 * Syntax:
 *   [Jti,x] = register_Jtinvw(J,detJ,w)
 * where
 *   J is a n-by-n-by-[sz] array
 *   detJ is an array of size [sz], each element detJ(i) containing
 *      the determinant of J(:,:,i)
 *   w a n-by-[sz] array of right-hand-side vectors
 * and
 *   Jti is a n-by-n-by-[sz] array, where Jti(:,:,i) = inv(J(:,:,i)');
 *   x is n-by-[sz] array, where x(:,i) is a solution to
 *             J(:,:,i)' * x(:,i) = w(:,i)
 *      (note the transpose of J!)
 *
 * Supplying detJ speeds things up, and we've already pre-computed it.
 *
 * Copyright 2006 by Timothy E. Holy
 */

void Jtinvwwork(const float *J,const float *Jend,int matrix_dim,const float *detJ,const float *w,float *Jti,float *x);
int validate_dimensions(const int *,int,const int *,int);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *J,*detJ,*w;
  float *Jti,*x;
  int n_dims_J,n_dims_detJ,n_dims_w,matrix_dim;
  const int *sz_J,*sz_detJ,*sz_w;
  const mxArray *curarg;

  if (nrhs != 3)
    mexErrMsgTxt("Jtinvw: requires three inputs");
  if (nlhs != 1)
    mexErrMsgTxt("Jtinvw: requires one output");

  // Parse the inputs
  // J:
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_Jtinvw: J must be a real single-precision array");
  J = (float *) mxGetData(curarg);
  n_dims_J = mxGetNumberOfDimensions(curarg);
  if (n_dims_J < 2)
    mexErrMsgTxt("register_Jtinvw: J must have at least two dimensions");
  sz_J = mxGetDimensions(curarg);
  // Sanity check on dimensionality
  matrix_dim = sz_J[0];
  if (matrix_dim != sz_J[1])
    mexErrMsgTxt("register_Jtinvw: J must be a n-by-n-by-[sz] array");
  if (matrix_dim > 3)
    mexErrMsgTxt("register_Jtinvw: for now, handles up to 3-by-3 matrices only");

  // detJ:
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_Jtinvw: J must be a real single-precision array");
  detJ = (float *) mxGetData(curarg);
  n_dims_detJ = mxGetNumberOfDimensions(curarg);
  sz_detJ = mxGetDimensions(curarg);
  if (!validate_dimensions(sz_J+2,n_dims_J-2,sz_detJ,n_dims_detJ))
    mexErrMsgTxt("register_Jtinvw: dimensions of J and detJ do not match");
    
  // w:
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_Jtinvw: w must be a real single-precision array");
  w = (float *) mxGetData(curarg);
  n_dims_w = mxGetNumberOfDimensions(curarg);
  sz_w = mxGetDimensions(curarg);
  if (!validate_dimensions(sz_J+1,n_dims_J-1,sz_w,n_dims_w))
    mexErrMsgTxt("register_Jtinvw: dimensions of J and w do not match");
    

  // Set up the outputs
  //plhs[0] = mxCreateNumericArray(n_dims_J,sz_J,mxSINGLE_CLASS,mxREAL);
  //plhs[1] = mxCreateNumericArray(n_dims_w,sz_w,mxSINGLE_CLASS,mxREAL);
  //Jti = (float *) mxGetData(plhs[0]);
  //x = (float *) mxGetData(plhs[1]);
  plhs[0] = mxCreateNumericArray(n_dims_w,sz_w,mxSINGLE_CLASS,mxREAL);
  x = (float *) mxGetData(plhs[0]);
  Jti = NULL;

  // Do the actual work
  Jtinvwwork(J,J+mxGetNumberOfElements(prhs[0]),matrix_dim,detJ,w,Jti,x);

  return;
}


/*
 * Do the actual computation.
 * It uses explicit formulas for inverse-transpose in dimensions 1, 2, and 3. 
 * Presumably the overhead for LU-decomposition would not be justifiable.
 */
void Jtinvwwork(const float *J,const float *Jend,int matrix_dim,const float *detJ,const float *w,float *Jti,float *x)
{
  int pixelSkip,Jti_pixelskip;
  float detJtmp;
  float Jtitmp[9];

  pixelSkip = matrix_dim*matrix_dim;
  Jti_pixelskip = pixelSkip;
  if (Jti == NULL) {
    Jti_pixelskip = 0;
    Jti = &(Jtitmp[0]);
  }
  if (matrix_dim == 1)
    for (; J < Jend; J++,w++,Jti+=Jti_pixelskip,x++) {
      *Jti = 1.0 / *J;
      *x = *w * *Jti;
    }
  else if (matrix_dim == 2)
    for (; J < Jend; J += pixelSkip, detJ++, w += matrix_dim, Jti += Jti_pixelskip, x += matrix_dim) {
      detJtmp = *detJ;
      Jti[0] = J[3] / detJtmp;
      Jti[1] = -J[2] / detJtmp;
      Jti[2] = -J[1] / detJtmp;
      Jti[3] = J[0] / detJtmp;
      x[0] = Jti[0]*w[0] + Jti[2]*w[1];
      x[1] = Jti[1]*w[0] + Jti[3]*w[1];
    }
  else if (matrix_dim == 3)
    for (; J < Jend; J += pixelSkip, detJ++, w += matrix_dim, Jti += Jti_pixelskip, x += matrix_dim) {
      detJtmp = *detJ;
      Jti[0] = (J[4]*J[8]-J[7]*J[5]) / detJtmp;
      Jti[3] = -(J[1]*J[8]-J[7]*J[2]) / detJtmp;
      Jti[6] = (J[1]*J[5]-J[4]*J[2]) / detJtmp;
      Jti[1] = -(J[3]*J[8]-J[6]*J[5]) / detJtmp;
      Jti[4] = (J[0]*J[8]-J[6]*J[2]) / detJtmp;
      Jti[7] = -(J[0]*J[5]-J[3]*J[2]) / detJtmp;
      Jti[2] = (J[3]*J[7]-J[6]*J[4]) / detJtmp;
      Jti[5] = -(J[0]*J[7]-J[6]*J[1]) / detJtmp;
      Jti[8] = (J[0]*J[4]-J[3]*J[1]) / detJtmp;
      x[0] = Jti[0]*w[0] + Jti[3]*w[1] + Jti[6]*w[2];
      x[1] = Jti[1]*w[0] + Jti[4]*w[1] + Jti[7]*w[2];
      x[2] = Jti[2]*w[0] + Jti[5]*w[1] + Jti[8]*w[2];
    }

  return;
}

int validate_dimensions(const int *d1,int l1,const int *d2, int l2)
{
  int i;

  if (l1 != l2)
    return 0;
  for (i = 0; i < l1; i++)
    if (d1[i] != d2[i])
      return 0;

  return 1;
}
