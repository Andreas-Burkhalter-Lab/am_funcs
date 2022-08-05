/*=================================================================
  Fill in the entries of an interpolation matrix

  Copyright 2006 by Timothy E. Holy
 *=================================================================*/

#include "mex.h"
//#include <vector>

//using namespace std;


//---------------------------------------------------------------------

void im_work(int n_dims,double **coords_in,double **coords_out,int *length_in,int *length_out,double *index_in,double *index_out,double *value);

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int n_dims,n_nbrs,n_pts_out;
  double **coords_in,**coords_out;
  int *length_in,*length_out;
  double *index_in,*index_out,*value;

  int carg,i;

  // Check for proper number of input and output arguments
  n_dims = nrhs/2;
  if (2*n_dims != nrhs) {
    mexErrMsgTxt("Input dimensionality doesn't match output dimensionality.");
  } 
  if(nlhs != 3)
    mexErrMsgTxt("Three output arguments required.");

  // Allocate space for storing data
  coords_in = (double **) mxMalloc(n_dims * sizeof(double*));
  coords_out = (double **) mxMalloc(n_dims * sizeof(double*));
  length_in = (int *) mxMalloc(n_dims * sizeof(int));
  length_out = (int *) mxMalloc(n_dims * sizeof(int));

  // Input coordinates
  for (carg = 0; carg < n_dims; carg++) {
    if( !mxIsDouble(prhs[carg]) ) mexErrMsgTxt("Input must be double.");
    if( mxIsComplex(prhs[carg]) ) mexErrMsgTxt("Input cannot be complex.");
    length_in[carg] = mxGetNumberOfElements(prhs[carg]);
    coords_in[carg] = mxGetPr(prhs[carg]);
  }
  // Output coordinates
  for (carg = 0; carg < n_dims; carg++) {
    if( !mxIsDouble(prhs[carg+n_dims]) ) mexErrMsgTxt("Input must be double.");
    if( mxIsComplex(prhs[carg+n_dims]) ) mexErrMsgTxt("Input cannot be complex.");
    length_out[carg] = mxGetNumberOfElements(prhs[carg+n_dims]);
    coords_out[carg] = mxGetPr(prhs[carg+n_dims]);
  }
  
  // Prepare space for output
  n_pts_out = 1;
  for (i = 0; i < n_dims; i++) {
    n_pts_out *= length_out[i];
  }
  n_nbrs = (1 << n_dims);
  plhs[0] = mxCreateDoubleMatrix(1,n_pts_out*n_nbrs,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,n_pts_out*n_nbrs,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,n_pts_out*n_nbrs,mxREAL);
  index_out = mxGetPr(plhs[0]);
  index_in = mxGetPr(plhs[1]);
  value = mxGetPr(plhs[2]);

  im_work(n_dims,coords_in,coords_out,length_in,length_out,index_in,index_out,value);

  mxFree(coords_in);
  mxFree(coords_out);
  mxFree(length_in);
  mxFree(length_out);
}

void im_work(int n_dims,double **coords_in,double **coords_out,int *length_in,int *length_out,double *index_in,double *index_out,double *value)
{
  double **f;
  int **ilr, *cumsize_in, *cumsize_out, *bitmask, *tcoords_out;
  int i,j,k,n_nbrs,ioffset_in,tindex,itmp,inbr;
  double fac;
  
  // Allocate all temporary storage
  ilr = (int **) mxMalloc(n_dims*sizeof(int *));
  f = (double **) mxMalloc(n_dims*sizeof(double *));
  for (i = 0; i < n_dims; i++) {
    ilr[i] = (int *) mxMalloc(2*length_out[i]*sizeof(int));
    f[i] = (double *) mxMalloc(length_out[i]*sizeof(double));
  }
  cumsize_in = (int *) mxMalloc((n_dims+1)*sizeof(int));
  cumsize_out = (int *) mxMalloc((n_dims+1)*sizeof(int));
  bitmask = (int *) mxMalloc(n_dims*sizeof(int));
  tcoords_out = (int *) mxMalloc(n_dims*sizeof(int));
  n_nbrs = (1 << n_dims);

  // Initialize the straightforward temporary storage
  cumsize_in[0] = cumsize_out[0] = 1;
  for (i = 1; i <= n_dims; i++) {
    cumsize_in[i] = cumsize_in[i-1]*length_in[i-1];
    cumsize_out[i] = cumsize_out[i-1]*length_out[i-1];
  }
  for (i = 0; i < n_dims; i++)
    bitmask[i] = (1 << i);
  n_nbrs = (1 << n_dims);

  // Loop over the dimensions
  for (i = 0; i < n_dims; i++) {
    // For each output coordinate, find the input index to the left
    tindex = -1;
    for (j = 0; j < length_out[i]; j++) {
      while (tindex+1 < length_in[i] && coords_in[i][tindex+1] <= coords_out[i][j])
	tindex++;
      ilr[i][2*j] = tindex;
    }
    // For each output coordinate, find the input index to the right
    tindex = length_in[i];
    for (j = length_out[i]-1; j > -1; j--) {
      while (tindex > 0 && coords_in[i][tindex-1] >= coords_out[i][j])
	tindex--;
      ilr[i][2*j+1] = tindex;
    }
    // Compute the linear factors
    for (j = 0; j < length_out[i]; j++) {
      if (ilr[i][2*j] != ilr[i][2*j+1])
	// The left & right points are different
	f[i][j] = (coords_out[i][j] - coords_in[i][ilr[i][2*j]]) / (coords_in[i][ilr[i][2*j+1]] - coords_in[i][ilr[i][2*j]]);
      else
	// The left & right are the same point
	f[i][j] = 1;
    }
  }

  // Now start filling output values
  for (i = 0; i < cumsize_out[n_dims]; i++) {
    // Convert output index to subscripts
    itmp = i;
    for (j = n_dims-1; j >= 0; j--) {
      tcoords_out[j] = itmp/cumsize_out[j];
      itmp -= tcoords_out[j] * cumsize_out[j];
    }
    // Loop over neighbors
    for (j = 0; j < n_nbrs; j++, index_in++, index_out++, value++) {
      ioffset_in = 0;
      fac = 1;
      for (k = 0; k < n_dims; k++) {
	inbr = ilr[k][2*tcoords_out[k] + ((j & bitmask[k])>0)]; // in-index of jth nbr
	ioffset_in += inbr*cumsize_in[k];
	if (j & bitmask[k])
	  fac *= f[k][tcoords_out[k]];
	else
	  fac *= 1.0 - f[k][tcoords_out[k]];
      }
      *index_in = ioffset_in+1;   // +1 for Matlab offset
      *index_out = i+1;           // ditto
      *value = fac;
    }
  }

  // Free temporary storage
  for (i = 0; i < n_dims; i++) {
    mxFree(ilr[i]);
    mxFree(f[i]);
  }
  mxFree(ilr);
  mxFree(f);
  mxFree(cumsize_in);
  mxFree(cumsize_out);
  mxFree(bitmask);
  mxFree(tcoords_out);
}
