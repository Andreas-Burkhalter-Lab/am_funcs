#include "mex.h"

// Matlab calling syntax: tridiag_inv(ldiag,diag,udiag,rhs)

template<class dataType>
void tridiag_inv_work(dataType *ldiag,dataType *diag,dataType *udiag,dataType *rhs,dataType *result,int M,int N,int matrix_inc_flag)
{
  dataType *m, k;
  int i,j;

  m = (dataType*) mxMalloc(M*sizeof(dataType));
  for (i = 0; i < N; i++, ldiag += matrix_inc_flag*(M-1), diag += matrix_inc_flag*M, udiag += matrix_inc_flag*(M-1), rhs += M, result += M) {
    k = *diag;
    *result = *rhs/k;
    for (j = 1; j < M; j++) {
      m[j] = udiag[j-1]/k;
      k = diag[j] - ldiag[j-1]*m[j];
      if (k == 0)
	mexErrMsgTxt("Tridiagonal inversion failed");
      result[j] = (rhs[j]-ldiag[j-1]*result[j-1])/k;
    }
    for (j = M-2; j >= 0; j--)
      result[j] -= m[j+1]*result[j+1];
  }
  mxFree(m);
}

template<class dataType>
void mexWrapper(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[],mxClassID matDataType)
{
  // diag is the diagonal, ldiag is the lower diagonal, and udiag is the upper diagonal
  // rhs is the rhs
  dataType *ldiag, *diag, *udiag, *rhs, *result;
  int M,N,n,matrix_inc_flag;
  const mxArray *curarg,*rhsarg;

  // Parse the inputs
  // Start with the rhs to get sizes
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("tridiag_inv: rhs must be a real matrix");
  M = mxGetM(curarg);
  N = mxGetNumberOfElements(curarg)/M;
  rhs = (dataType *) mxGetData(curarg);
  rhsarg = curarg;
  
  // ldiag
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("tridiag_inv: ldiag must be a real matrix");
  ldiag = (dataType *) mxGetData(curarg);
  if (mxGetNumberOfElements(curarg) != M-1 && mxGetM(curarg) != M-1)
    mexErrMsgTxt("tridiag_inv: ldiag must be a vector of length M-1 or have columns of size M-1");
  n = mxGetNumberOfElements(curarg)/(M-1);
  if (n == 1)
    matrix_inc_flag = 0;
  else if (n == N)
    matrix_inc_flag = 1;
  else
    mexErrMsgTxt("tridiag_inv: ldiag must be a vector or have its remaining dimensions be the same size as rhs");

  // diag
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("tridiag_inv: diag must be a real matrix");
  if (mxGetNumberOfElements(curarg) != M || mxGetNumberOfElements(curarg)/M != n) {
    mexPrintf("M %d, n %d, numel %d\n",M,n,mxGetNumberOfElements(curarg));
    mexErrMsgTxt("tridiag_inv: diag must be a vector of length M or have columns of size M, and its remaining dimensions of the same size as ldiag");
  }
  diag = (dataType *) mxGetData(curarg);

  // udiag
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("tridiag_inv: udiag must be a real matrix");
  if (mxGetNumberOfElements(curarg) != M-1 || mxGetNumberOfElements(curarg)/(M-1) != n)
    mexErrMsgTxt("tridiag_inv: udiag must be a vector of length M-1 or have columns of size M-1, and its remaining dimensions of the same size as ldiag");
  udiag = (dataType *) mxGetData(curarg);

  // Set up the output
  plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(rhsarg),mxGetDimensions(rhsarg),matDataType,mxREAL);
  result = (dataType *) mxGetData(plhs[0]);

  // Do the actual work
  tridiag_inv_work(ldiag,diag,udiag,rhs,result,M,N,matrix_inc_flag);
}


void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mxClassID matDataType;
  int i;

  if (nrhs != 4)
    mexErrMsgTxt("tridiag_inv: requires four inputs");

  // Check the data types for consistency
  matDataType = mxGetClassID(prhs[0]);
  for (i = 1; i < 4; i++)
    if (mxGetClassID(prhs[i]) != matDataType)
      mexErrMsgTxt("tridiag_inv: all inputs must have the same data type");

  if (matDataType == mxSINGLE_CLASS)
    mexWrapper<float>(nlhs,plhs,nrhs,prhs,matDataType);
  else if (matDataType == mxDOUBLE_CLASS)
    mexWrapper<double>(nlhs,plhs,nrhs,prhs,matDataType);
  else
    mexErrMsgTxt("tridiag_inv: all inputs must be either single or double");
}


