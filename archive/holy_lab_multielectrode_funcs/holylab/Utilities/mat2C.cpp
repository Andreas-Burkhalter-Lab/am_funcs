// Some convenience functions for working with matlab structures. Note that these routines do a fair amount of error checking, but are not arbitrarily robust; basically, they guard against memory errors, but do not guard against float/int32 conversions, etc. The idea is to not get in the way of performance by having to cast each value.

template <class T>
    void matScalarField2C(const mxArray *mxPtr,const char *name,T *v)
{ 
  mxArray *fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr == NULL)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: missing required field '%s'",mexFunctionName(),name);
  if (mxGetNumberOfElements(fieldPtr) != 1)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: expect field '%s' to be a scalar",mexFunctionName(),name);
  *v = (T) mxGetScalar(fieldPtr);
}

template <class T>
    void matOptionalScalarField2C(const mxArray *mxPtr,const char *name,T *v)
{
  const mxArray *fieldPtr;

  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr != NULL) {
    if (mxGetNumberOfElements(fieldPtr) != 1)
      mexErrMsgIdAndTxt("custommex:field parsing error","%s: expect field '%s' to be a scalar",mexFunctionName(),name);
    *v = (T) mxGetScalar(fieldPtr);
  }
}

void C2matScalarField(mxArray *mxPtr,const char *name,double v)
{
  mxArray *fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr == NULL)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: required field '%s' not found",mexFunctionName(),name);
  if (!mxIsDouble(fieldPtr))
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' must be of type double",mexFunctionName(),name);
  double *p = mxGetPr(fieldPtr);
  *p = v;
}

// You cannot "cast" with the following function, it simply does a memcpy
template <class T>
    int matField2Carray(const mxArray *mxPtr,const char *name,T *v,int nmax,int nmin)
{
  const mxArray *fieldPtr;

  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr == NULL)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: required field '%s' not found",mexFunctionName(),name);
  const int n_bytes = mxGetElementSize(fieldPtr);
  if (n_bytes != sizeof(T))
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' is being converted to a type with the wrong number of bytes (%d instead of its native %d)",mexFunctionName(),name,sizeof(T),n_bytes);
  int n = mxGetNumberOfElements(fieldPtr);
  if (n > nmax)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' has %d elements, which is larger than the maximum number %d",mexFunctionName(),name,n,nmax);
  if (n < nmin)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' has %d elements, which is smaller than the minimum number %d",mexFunctionName(),name,n,nmax);
  T *data = (T*) mxGetData(fieldPtr);
  memcpy(v,data,n*sizeof(T));
  return n;
}

// Your array already must be allocated, of the right type, and of the proper size/shape. This simply does a memcpy
template<class T>
    void Carray2matField(mxArray *mxPtr,const char *name,const T *v,int n)
{
  const mxArray *fieldPtr;

  fieldPtr = mxGetField(mxPtr,0,name);
  if (fieldPtr == NULL)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: required field '%s' not found",mexFunctionName(),name);
  const int n_bytes = mxGetElementSize(fieldPtr);
  if (n_bytes != sizeof(T))
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' is being converted to a type with the wrong number of bytes (%d instead of its native %d)",mexFunctionName(),name,sizeof(T),n_bytes);
  int nelem = mxGetNumberOfElements(fieldPtr);
  if (nelem != n)
    mexErrMsgIdAndTxt("custommex:field parsing error","%s: field '%s' has %d elements, which is not equal to the desired number (%d)",mexFunctionName(),name,nelem,n);
  T *data = (T*) mxGetData(fieldPtr);
  memcpy(data,v,n*sizeof(T));
}
