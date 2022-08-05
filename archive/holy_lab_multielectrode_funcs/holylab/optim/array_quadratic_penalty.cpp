//#include <iostream>
#include "Eigen/Dense"
#include "imiterators.cxx"
#include "mex.h"
//#include "mat.h"
#include "MatlabIO.h"

using namespace std;
using namespace Eigen;

// Syntax:
//   [v,g] = array_quadratic_penalty(A,vertex,Q,boundary_conditions)
// where
//   A is a d-dimensional array;
//   vertex is a 2^d-by-d matrix of 0s and 1s, giving the displacement
//     of each "corner" of a unit cell (this matrix defines how Q is
//     interpreted);
//   Q is a 2^d-by-2^d matrix, specifying the penalty on a single unit cell:
//        Pcell = (1/2) * a^T * Q * a,
//     where a would be a vector of length 2^d representing the values
//     of A at the vertices of a single unit cell, in the order
//     specified by the variable "vertex";
//   boundary_conditions (optional) is a string, one of 'none'
//     (default) or 'circular', specifying how boundaries are supposed
//     to be treated.  'none' means that the penalty is applied only
//     to cells that lie fully within the bounds of the array; 'zero'
//     (not yet implemented) assumes that the function values are zero
//     beyond the edge of the array; 'circular' uses periodic boundary
//     conditions.
// On output,
//   v is the value of the penalty, summed over all unit cells (a scalar);
//   g is the gradient of v with respect to the coordinates of A, and
//     therefore is an array of the size of A.
//
// For example, in two dimensions, if we wanted to use a penalty function
//     Pcell = (1/2) (grad A)^2,
// then
//   vertex = [0 0;
//             0 1;
//             1 0;
//             1 1];
// which indicates a unit cell vertex numbering scheme like this:
//
//        3     4
//
//        1     2
//
// (In this coordinate system, we'll call the second coordinate "x",
// which increases to the right; y increases upward.)
// The x-gradient is
//       (a(2)+a(4)-a(1)-a(3))/2
// (the average value on the two right vertices, minus the average
// value on the two left vertices), and the y-gradient is
//       (a(3)+a(4)-a(1)-a(2))/2.
// Let gx = [-1/2, 1/2, -1/2, 1/2]^T, gy = [-1/2,-1/2,1/2,1/2]^T, in
// terms of which the gradients are expressed as a dot product.  We
// have
//       Q = gx gx^T + gy gy^T,
// and hence
//            [ 1/2   0    0  -1/2 ] 
//      Q =   [  0   1/2 -1/2   0  ]
//            [  0  -1/2  1/2   0  ]
//            [-1/2   0    0   1/2 ]

// Copyright 2011 by Timothy E. Holy

template <typename DataType>
struct TempSpace {
  TempSpace(int n_dims) {
    int n_vertices = 1 << n_dims;
    offset.resize(n_vertices);
    offset2.resize(n_vertices);
    cellv.resize(n_vertices);
    Qcellv.resize(n_vertices);
  }

  VectorXi offset,offset2;
  Matrix<DataType,Dynamic,1> cellv,Qcellv;
};

enum BCType {
  BCNONE,
  BCPER
};

template <typename Derived>
void aqp_offset(VectorXi &offset,const MatrixBase<Derived> &vertex, const pixIterator &pixI)
{
  for (int vertexIter = 0; vertexIter < offset.size(); vertexIter++) {
    int thisOffset = 0;
    for (int dimIter = 0; dimIter < pixI.nDims(); dimIter++) {
      int thisCoord = pixI.coord(dimIter) + vertex(vertexIter,dimIter);
      thisCoord = thisCoord % pixI.dimSize(dimIter);  // periodic BCs
      thisOffset += pixI.dimSkip(dimIter) * thisCoord;
    }
    offset[vertexIter] = thisOffset;
  }
  
  //cout << "Offset for " << pixI.coord(0) << ' ' << pixI.coord(1) << ": " << offset.transpose() << endl;
}

template <typename Derived>
void aqp_work(typename Derived::Scalar &v,typename Derived::Scalar *g,const typename Derived::Scalar *A,pixIterator &pixI,const MatrixBase<Derived>& vertex,const MatrixBase<Derived>& Q,BCType mode,TempSpace<typename Derived::Scalar>& tmp)
{
  int pixIter;

  // Compute the default memory offset vector
  aqp_offset(tmp.offset,vertex,pixI);
  
  v = 0;
  for (; !pixI.at_end(); pixI++) {
    if (!pixI.on_right_edge()) {
      for (pixIter = 0; pixIter < tmp.offset.size(); pixIter++)
	tmp.cellv[pixIter] = A[pixI+tmp.offset[pixIter]];
    } else {
      if (mode == BCNONE)
	continue;
      else if (mode == BCPER) {
	// Compute the offset appropriate to this particular unit cell
	aqp_offset(tmp.offset2,vertex,pixI);
	for (pixIter = 0; pixIter < tmp.offset.size(); pixIter++)
	  tmp.cellv[pixIter] = A[tmp.offset2[pixIter]];
      }
    }
    tmp.Qcellv = Q*tmp.cellv;
    v += tmp.Qcellv.dot(tmp.cellv);
    if (g != NULL) {
      if (!pixI.on_right_edge())
	for (pixIter = 0; pixIter < tmp.offset.size(); pixIter++)
	  g[pixI+tmp.offset[pixIter]] += tmp.Qcellv[pixIter];
      else
	for (pixIter = 0; pixIter < tmp.offset.size(); pixIter++)
	  g[tmp.offset2[pixIter]] += tmp.Qcellv[pixIter];
    }
  }
  v = v/2;
}

template <typename DataType>
void mexWrapper(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
  typedef typename MatlabIO::MapMx<DataType>::MapMatrixType MatrixType;
  const DataType *A;
  int n_dims,n_dims1,two_dims;
  const int *sz;
  char modechar;
  const int sz_scalar[2] = {1, 1};
  const mxArray *curarg;
  int dimIter;

  // Parse the arguments
  // A
  curarg = prhs[0];
  A = (DataType*) mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  sz = mxGetDimensions(curarg);
  for (n_dims1 = 0, dimIter = 0; dimIter < n_dims; dimIter++)
    n_dims1 += (sz[dimIter] > 1);
  two_dims = 1<<n_dims1;  // 2^n_dims
  // vertex
  curarg = prhs[1];
  if (mxGetClassID(curarg) != MatlabTraits<DataType>::ClassID())
    mexErrMsgTxt("ClassID of vertex does not match that of A");
  if (mxGetM(curarg) != two_dims || mxGetN(curarg) != n_dims1)
    mexErrMsgTxt("Size of vertex is inconsistent with A");
  MatrixType vertex(NULL,0,0);
  MatlabIO::mx2eigen(curarg,vertex);
  // Q
  curarg = prhs[2];
  if (mxGetClassID(curarg) != MatlabTraits<DataType>::ClassID())
    mexErrMsgTxt("ClassID of Q does not match that of A");
  if (mxGetM(curarg) != two_dims || mxGetN(curarg) != two_dims)
    mexErrMsgTxt("Size of Q is inconsistent with A");
  MatrixType Q(NULL,0,0);
  MatlabIO::mx2eigen(curarg,Q);
  // boundary_conditions
  BCType mode = BCNONE;
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsChar(curarg))
      mexErrMsgTxt("boundary_conditions must be a character array");
    modechar = *((char *) mxGetData(curarg));  // get first character
    if (modechar == 'p' || modechar == 'P')
      mode = BCPER;
    else if (modechar == 'n' || modechar == 'N')
      ;
    else
      mexErrMsgTxt("Mode not recognized");
  }

  // Allocate space for the output
  plhs[0] = MatlabTraits<DataType>::mxAllocator(2,sz_scalar);
  DataType *vp = (DataType*) mxGetData(plhs[0]);
  DataType *g = NULL;
  if (nlhs > 1) {
    plhs[1] = MatlabTraits<DataType>::mxAllocator(n_dims,sz);
    g = (DataType*) mxGetData(plhs[1]);
  }

  DataType v;
  TempSpace<DataType> tmp(n_dims1);
  pixIterator pixI(sz,n_dims);

  aqp_work(v,g,A,pixI,vertex,Q,mode,tmp);

  *vp = v;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  bool isdouble;
  const mxArray *curarg;

  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("array_quadratic_penalty: requires three or four inputs");
  if (nlhs < 1 || nlhs > 2)
    mexErrMsgTxt("array_quadratic_penalty: requires one or two outputs");

  // Parse the input to determine data type
  // A
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("array_quadratic_penalty: A must be a real array");
  if (mxIsSingle(curarg))
    isdouble = false;
  else if (mxIsDouble(curarg))
    isdouble = true;
  else
    mexErrMsgTxt("array_quadratic_penalty: A must be single or double");
  if (mxIsEmpty(curarg))
    mexErrMsgTxt("array_quadratic_penalty: does not work on empty arrays");

  if (isdouble)
    mexWrapper<double>(nlhs,plhs,nrhs,prhs);
  else
    mexWrapper<float>(nlhs,plhs,nrhs,prhs);
}

/*
int main(int argc,char *argv[])
{
  if (argc < 3)
    exit(1);

  //typedef Matrix<double,Dynamic,Dynamic> MatrixType;

  MatlabIO::Load loader(argv[1]);
  MatlabIO::Save saver(argv[2]);
  MatlabIO::MapMxD A(NULL,0,0);
  loader.load("A",A);
  MatlabIO::MapMxD Q(NULL,0,0);
  loader.load("Q",Q);
  MatlabIO::MapMxD vertex(NULL,0,0);
  loader.load("vertex",vertex);

  double v;
  Matrix<double,Dynamic,Dynamic> g(A.rows(),A.cols());

  int sz[2];
  sz[0] = A.rows();
  sz[1] = A.cols();
  pixIterator pixI(sz,2);
  TempSpace<double> tmp(2);

  aqp_work(v,&g(0,0),&A(0,0),pixI,vertex,Q,BCNONE,tmp);

  saver.saveScalar(v,"v");
  saver.save(g,"g");
}
*/
