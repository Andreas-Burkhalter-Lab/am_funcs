#include "image_utilities.h"

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

int skip_unity_dimensions(const int *din,int n,int *dout,int noutmax)
{
  int i1,i2;

  i1 = 0;
  i2 = 0;
  while (i1 < n && i2 < noutmax) {
    dout[i2] = din[i1];
    if (din[i1] > 1)
      i2++;
    i1++;
  }
  if (i1 < n)
    return -1;  // Error, too many dimensions

  return i2;  // the number of dimensions
}

int skip_unity_dimensions_index(const int *din,int n,int *indxout,int noutmax)
{
  int i1,i2;

  i1 = 0;
  i2 = 0;
  while (i1 < n && i2 < noutmax) {
    indxout[i2] = i1;
    if (din[i1] > 1)
      i2++;
    i1++;
  }
  if (i1 < n)
    return -1;  // Error, too many dimensions

  return i2;
}

long calc_pixel_skip(const int *sz,int n_dims,long *pixel_skip)
{
  int dimIndex;
  long ret;

  // Figure out the memory spacing between adjacent pixels in all directions
  pixel_skip[0] = 1;
  for (dimIndex = 1; dimIndex < n_dims; dimIndex++)
    pixel_skip[dimIndex] = pixel_skip[dimIndex-1] * sz[dimIndex-1];
  ret = pixel_skip[n_dims-1]*sz[n_dims-1];  // the total size
  // For any dimensions with size=1, set pixel_skip to 0
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (sz[dimIndex] == 1)
      pixel_skip[dimIndex] = 0;

  return ret;
}

int isVector(const int *sz,const int n_dims)
{
  int n_nonunity,dimIndex;

  n_nonunity = 0;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (sz[dimIndex] > 1)
      n_nonunity++;
  return (n_nonunity < 2);
}
