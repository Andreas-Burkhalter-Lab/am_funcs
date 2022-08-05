#include <cmath>
#include <limits>

// qinterp: class for quadratic interpolation
// Usage:
//   Qinterp::qinterp<dataType,n_dims> q(sz)
// This is an example of the basic declaration, where dataType is a
// floating-precision type (float, double, or long double) and n_dims
// is 1, 2, or 3. Values over the edge of the array are given a value
// of NaN. sz is an array of integers, giving the dimensions of the
// n_dims-dimensional array on which to perform interpolation.
//
//   Qinterp::qinterp<dataType,n_dims,Qinterp::reflect,1> q(sz)
// A more customized declaration.  The third argument specifies how the
// boundaries are to be treated; options are nocheck (which skips all
// checks and is therefore fastest), NaN (values outside the array are
// NaN), reflect (using reflecting boundary conditions), and periodic
// (using periodic boundary conditions).
// The last template argument contains the coordinate offset; specify
// 1 if your coordinates are unit-offset (like in Matlab) or 0 if
// zero-offset.
//
// To perform interpolation:
//   q.value(x,A,&val)
//   q.valgrad(x,A,&val,gP)
// These two calculate either the value (first line) or value and
// gradient (second line) at a single point x, for a single array A.
// gP points to an array with at least n_dims elements.
//
//   q.value(x,3,AP,valP)
//   q.valgrad(x,3,AP,valP,gP)
// These perform interpolation for 3 arrays AP[0], AP[1], and AP[2] at
// a single point x. gP points to an array of at least 3*n_dims
// elements.  The order of values is
//   AP0grad0 AP1grad0 AP2grad0 AP0grad1 AP1grad1 ...
// where gradi means the derivative with respect to the ith coordinate.

// Copyright 2011 by Timothy E. Holy

namespace Qinterp {

#define DIMENSIONALITY_NOT_IMPLEMENTED

enum boundary_condition {nocheck, NaN, reflect, periodic};

// Compile-time computation of powers of 3
template <int n>
struct pow3
{
  enum { value = 3*pow3<n-1>::value };
};
template <>
struct pow3<0>
{
  enum { value = 1 };
};

//
// Quadratic interpolation class definition
//
template <class dataType,int n_dims,boundary_condition bc = NaN,int coord_offset = 0>
class qinterp
{
  int sz[n_dims];    // size of the array
  size_t cumsz[n_dims]; // memory offset for increments along coordinates
  int xint[n_dims];  // "base point", rounded coordinate values
  dataType coef1d[3*n_dims];   // coefficients for value along each dim.
  dataType coefg1d[3*n_dims];  // coefficients for gradient along each dim.
  dataType coefprod[pow3<n_dims>::value];  // product of coefficients
  size_t memoffset[pow3<n_dims>::value];   // memory offsets of neighbors (typical)
  size_t thisoffset[pow3<n_dims>::value];  // offsets of neighbors (this point), used only for certain boundary conditions

  // Calculate the 1-d coefficients
  void set_coef1d(const dataType x[]);  // for just value
  void set_coefg1d(const dataType x[]); // for value and gradient
  // Compute the products of coefficients across dimensions
  void set_coefprod();         // for value
  void set_coefgprod(int dim); // for gradient
  // Handling boundary conditions
  bool edge_test();
  void edge_apply();
public:
  // Customizing to a particular array size
  void set_array_size(const int szi[]);

  // Constructors
  qinterp() {;}
  qinterp(const int szi[]) {set_array_size(szi);}

  // Computation of interpolated value: a single array
  void value(const dataType x[],const dataType *A,dataType *val);
  // Computation of interpolated value: multiple arrays at the same point
  void value(const dataType x[],int n,const dataType *A[],dataType *val);

  // Computation of both value and the gradient: a single array
  void valgrad(const dataType x[],const dataType *A,dataType *val,dataType *grad);
  // Computation of both value and the gradient: multiple arrays
  void valgrad(const dataType x[],int n,const dataType *A[],dataType *val,dataType *grad);
};


// Use specialization to handle some aspects of the dimensionality; in
// some ways, it's easier than writing arbitrary-dimension code. But
// most importantly, it is guaranteed to run fast.
// (Note that by including the dimensionality as a template parameter,
// the test is compiled out, so there is no overhead)

// set_array_size
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::set_array_size(const int szi[])
{
  if (n_dims == 1) {
    sz[0] = szi[0];
    cumsz[0] = 1;
    memoffset[0] = -1;
    memoffset[1] = 0;
    memoffset[2] = 1;
  } else if (n_dims == 2) {
    sz[0] = szi[0];
    sz[1] = szi[1];
    cumsz[0] = 1;
    cumsz[1] = sz[0];
    int i = 0;
    for (int dimIndex1 = -1; dimIndex1 < 2; dimIndex1++)
      for (int dimIndex2 = -1; dimIndex2 < 2; dimIndex2++)
	memoffset[i++] = dimIndex1 + dimIndex2*cumsz[1];
  } else if (n_dims == 3) {
    sz[0] = szi[0];
    sz[1] = szi[1];
    sz[2] = szi[2];
    cumsz[0] = 1;
    cumsz[1] = sz[0];
    cumsz[2] = sz[1]*cumsz[1];
    int dimIndex[3];
    int i = 0;
    for (dimIndex[0] = -1; dimIndex[0] < 2; dimIndex[0]++)
      for (dimIndex[1] = -1; dimIndex[1] < 2; dimIndex[1]++)
	for (dimIndex[2] = -1; dimIndex[2] < 2; dimIndex[2]++) {
	  size_t tmpv = 0;
	  for (int j = 0; j < 3; j++)
	    tmpv += dimIndex[j]*cumsz[j];
	  memoffset[i++] = tmpv;
	}
  } else {
    // Arbitrary-dimension code
    int dimIndex,nbrIndex;
    size_t offset;
    bool carry;

    cumsz[0] = 1;
    xint[0] = -1;  // use xint as coordinate offset
    for (dimIndex = 1; dimIndex < n_dims; dimIndex++) {
      cumsz[dimIndex] = cumsz[dimIndex-1]*sz[dimIndex-1];
      xint[dimIndex] = -1;
    }
    xint[n_dims] = 1;  // sentinel value
    nbrIndex = 0;
    carry = false;
    while (!carry) {
      // Calculate the total memory offset at the current coordinate displacement
      offset = 0;
      for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
	offset = xint[dimIndex]*cumsz[dimIndex];
      memoffset[nbrIndex++] = offset;
      // Go on to the next neighbor
      xint[0]++;
      dimIndex = 0;
      carry = (xint[0] == 2);
      while (carry) {
	xint[dimIndex++] = -1;
	xint[dimIndex]++;
	carry = (xint[dimIndex] == 2);
      }
    }
  }
}

// set_coef1d
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::set_coef1d(const dataType x[])
{
  dataType xf,tmp;
  int i = 0;
  for (int dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    xint[dimIndex] = (int) floor(x[dimIndex]+0.5);
    xf = x[dimIndex] - xint[dimIndex];
    xint[dimIndex] = xint[dimIndex]-coord_offset;
    for (int j = -1; j < 2; j++) {
      if (j == 0)
	coef1d[i++] = 3.0/4-xf*xf;
      else {
	tmp = xf + dataType(j)/2;
	coef1d[i++] = tmp*tmp/2;
      }
    }
  }
}

// set_coefg1d
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::set_coefg1d(const dataType x[])
{
  dataType xf,tmp;
  int i = 0;
  for (int dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    xint[dimIndex] = (int) floor(x[dimIndex]+0.5);
    xf = x[dimIndex] - xint[dimIndex];
    xint[dimIndex] = xint[dimIndex]-coord_offset;
    for (int j = -1; j < 2; j++) {
      if (j == 0) {
	coef1d[i] = 3.0/4-xf*xf;
	coefg1d[i++] = -2*xf;
      }
      else {
	tmp = xf + dataType(j)/2;
	coef1d[i] = tmp*tmp/2;
	coefg1d[i++] = tmp;
      }
    }
  }
}

// edge_test
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
bool qinterp<dataType,n_dims,bc,coord_offset>::edge_test()
{
  // Check to see whether we will need to fiddle with the edges
  for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (xint[dimIndex] <= 0 || xint[dimIndex] >= sz[dimIndex]-1)
      return true;
  return false;
}


// value
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::value(const dataType x[],const dataType *A,dataType *val)
{
  dataType const *Av[1] = {A};

  value(x,1,Av,val);
}


template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::value(const dataType x[],int n,dataType const *A[],dataType *val)
{
  // When appropriate, check array bounds for an early exit
  if (bc == NaN) {
    for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
      if (x[dimIndex] < coord_offset+0.5 || x[dimIndex] >= coord_offset+sz[dimIndex]-1.5) {
	for (int outIndex = 0; outIndex < n; outIndex++)
	  val[outIndex] = std::numeric_limits<dataType>::quiet_NaN();
	return;
      }
  }
  // Calculate the interpolation coefficients
  set_coef1d(x);
  set_coefprod();
  // Calculate the value
  if ((bc == reflect || bc == periodic) && edge_test()) {
    edge_apply();
    for (int j = 0; j < n; j++) {
      dataType tmp = 0;
      for (int i = 0; i < pow3<n_dims>::value; i++)
	tmp += coefprod[i]*A[j][thisoffset[i]];
      val[j] = tmp;
    }
  } else {
    size_t xoffset = 0;
    for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
      xoffset += xint[dimIndex]*cumsz[dimIndex];
    for (int j = 0; j < n; j++) {
      const dataType *Atmp = A[j]+xoffset;
      dataType tmp = 0;
      for (int i = 0; i < pow3<n_dims>::value; i++)
	tmp += coefprod[i]*Atmp[memoffset[i]];
      val[j] = tmp;
    }
  }
}

// valgrad
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::valgrad(const dataType x[],const dataType *A,dataType *val,dataType *grad)
{
  dataType const *Av[1] = {A};

  valgrad(x,1,Av,val,grad);
}

template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::valgrad(const dataType x[],int n,const dataType *A[],dataType *val,dataType *grad)
{
  // When appropriate, check array bounds for an early exit
  if (bc == NaN) {
    for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
      if (x[dimIndex] < coord_offset+0.5 || x[dimIndex] >= coord_offset+sz[dimIndex]-1.5) {
	for (int outIndex = 0; outIndex < n; outIndex++)
	  val[outIndex] = std::numeric_limits<dataType>::quiet_NaN();
	return;
      }
  }
  // Calculate the interpolation coefficients
  set_coefg1d(x);
  // Calculate the value & gradient
  if ((bc == reflect || bc == periodic) && edge_test()) {
    edge_apply(); // case where we need to apply boundary conditions
    set_coefprod();
    for (int j = 0; j < n; j++) {
      dataType tmp = 0;
      for (int i = 0; i < pow3<n_dims>::value; i++)
	tmp += coefprod[i]*A[j][thisoffset[i]];
      val[j] = tmp;
    }
    int index = 0;
    for (int k = 0; k < n_dims; k++) {
      set_coefgprod(k);
      for (int j = 0; j < n; j++) {
	dataType tmp = 0;
	for (int i = 0; i < pow3<n_dims>::value; i++)
	  tmp += coefprod[i]*A[j][thisoffset[i]];
	grad[index++] = tmp;
      }
    }
  } else { // default case where point is in the interior
    size_t xoffset = 0;
    for (int dimIndex = 0; dimIndex < n_dims; dimIndex++)
      xoffset += xint[dimIndex]*cumsz[dimIndex];
    set_coefprod();
    for (int j = 0; j < n; j++) {
      const dataType *Atmp = A[j]+xoffset;
      dataType tmp = 0;
      for (int i = 0; i < pow3<n_dims>::value; i++)
	tmp += coefprod[i]*Atmp[memoffset[i]];
      val[j] = tmp;
    }
    int index = 0;
    for (int k = 0; k < n_dims; k++) {
      set_coefgprod(k);
      for (int j = 0; j < n; j++) {
	const dataType *Atmp = A[j]+xoffset;
	dataType tmp = 0;
	for (int i = 0; i < pow3<n_dims>::value; i++)
	  tmp += coefprod[i]*Atmp[memoffset[i]];
	grad[index++] = tmp;
      }
    }
  }
}


// set_coefprod
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::set_coefprod()
{
  if (n_dims == 1) {
    coefprod[0] = coef1d[0];
    coefprod[1] = coef1d[1];
    coefprod[2] = coef1d[2];
  } else if (n_dims == 2) {
    int dimIndex0,dimIndex1;
    int i = 0;
    for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
      for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	coefprod[i++] = coef1d[dimIndex0]*coef1d[dimIndex1+3];
  } else if (n_dims == 3) {
    int dimIndex0,dimIndex1,dimIndex2;
    int i = 0;
    for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
      for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	for (dimIndex2 = 0; dimIndex2 < 3; dimIndex2++)
	  coefprod[i++] = coef1d[dimIndex0]*coef1d[dimIndex1+3]*coef1d[dimIndex2+6];
  } else {
    DIMENSIONALITY_NOT_IMPLEMENTED;
  }
}

// set_coefgprod
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::set_coefgprod(int dim)
{
  if (n_dims == 1) {
    coefprod[0] = coefg1d[0];
    coefprod[1] = coefg1d[1];
    coefprod[2] = coefg1d[2];
  } else if (n_dims == 2) {
    int dimIndex0,dimIndex1;
    int i = 0;
    switch (dim) {
    case 0:
      for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
	for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	  coefprod[i++] = coefg1d[dimIndex0]*coef1d[dimIndex1+3];
      break;
    case 1:
      for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
	for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	  coefprod[i++] = coef1d[dimIndex0]*coefg1d[dimIndex1+3];
      break;
    }
  } else if (n_dims == 3) {
    int dimIndex0,dimIndex1,dimIndex2;
    int i = 0;
    switch (dim) {
    case 0:
      for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
	for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	  for (dimIndex2 = 0; dimIndex2 < 3; dimIndex2++)
	    coefprod[i++] = coefg1d[dimIndex0]*coef1d[dimIndex1+3]*coef1d[dimIndex2+6];
      break;
    case 1:
      for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
	for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	  for (dimIndex2 = 0; dimIndex2 < 3; dimIndex2++)
	    coefprod[i++] = coef1d[dimIndex0]*coefg1d[dimIndex1+3]*coef1d[dimIndex2+6];
      break;
    case 2:
      for (dimIndex0 = 0; dimIndex0 < 3; dimIndex0++)
	for (dimIndex1 = 0; dimIndex1 < 3; dimIndex1++)
	  for (dimIndex2 = 0; dimIndex2 < 3; dimIndex2++)
	    coefprod[i++] = coef1d[dimIndex0]*coef1d[dimIndex1+3]*coefg1d[dimIndex2+6];
      break;
    }
  } else {
    DIMENSIONALITY_NOT_IMPLEMENTED;
  }
}

// edge_apply
template <class dataType,int n_dims,boundary_condition bc,int coord_offset>
void qinterp<dataType,n_dims,bc,coord_offset>::edge_apply()
{
  if (n_dims == 1) {
    for (int nbrIndex = -1; nbrIndex < 2; nbrIndex++) {
      int thisxint = xint[0]+nbrIndex;
      if (bc == periodic) {
	thisxint = thisxint % sz[0];
	if (thisxint < 0)
	  thisxint += sz[0];
      }
      else if (bc == reflect) {
	int sz2 = 2*sz[0];
	if (thisxint < 0)
	  thisxint = -thisxint - 1;
	thisxint = thisxint % sz2;
	if (thisxint >= sz[0])
	  thisxint = sz2-thisxint-1;
      }
      thisoffset[nbrIndex+1] = thisxint;
    }
  } else if (n_dims == 2) {
    int nbrIndex[2];
    int i = 0;
    for (nbrIndex[0] = -1; nbrIndex[0] < 2; nbrIndex[0]++)
      for (nbrIndex[1] = -1; nbrIndex[1] < 2; nbrIndex[1]++) {
	size_t tmpv = 0;
	for (int j = 0; j < 2; j++) {
	  int thisxint = xint[j]+nbrIndex[j];
	  if (bc == periodic) {
	    thisxint = thisxint % sz[j];
	    if (thisxint < 0)
	      thisxint += sz[j];
	  }
	  else if (bc == reflect) {
	    int sz2 = 2*sz[j];
	    if (thisxint < 0)
	      thisxint = -thisxint - 1;
	    thisxint = thisxint % sz2;
	    if (thisxint >= sz[j])
	      thisxint = sz2-thisxint-1;
	  }
	  tmpv += thisxint*cumsz[j];
	}
	thisoffset[i++] = tmpv;
      }
  } else if (n_dims == 3) {
    int nbrIndex[3];
    int i = 0;
    for (nbrIndex[0] = -1; nbrIndex[0] < 2; nbrIndex[0]++)
      for (nbrIndex[1] = -1; nbrIndex[1] < 2; nbrIndex[1]++)
	for (nbrIndex[2] = -1; nbrIndex[2] < 2; nbrIndex[2]++) {
	  size_t tmpv = 0;
	  for (int j = 0; j < 3; j++) {
	    int thisxint = xint[j]+nbrIndex[j];
	    if (bc == periodic) {
	      thisxint = thisxint % sz[j];
	      if (thisxint < 0)
		thisxint += sz[j];
	    }
	    else if (bc == reflect) {
	      int sz2 = 2*sz[j];
	      if (thisxint < 0)
		thisxint = -thisxint - 1;
	      thisxint = thisxint % sz2;
	      if (thisxint >= sz[j])
		thisxint = sz2-thisxint-1;
	    }
	    tmpv += thisxint*cumsz[j];
	  }
	  thisoffset[i++] = tmpv;
	}

  } else {
    DIMENSIONALITY_NOT_IMPLEMENTED;
  }
}
}

//using Qinterp::qinterp;
