#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#else
#include "mex.h"
#endif

#include <math.h>   /* for floorf */
#include "imiterators.cxx"  /* for pixIterator */
#include "image_utilities.h"

/*
 * register_gradE: compute the gradient of the image difference error
 * with respect to coordinate deformation
 *
 * Syntax:
 *   grad = register_gradE(g,psi1,psi2,psig,sqrtdetJ)
 * where
 *   (pass empty sqrtdetJ, or leave input out, if you don't want the covariant theory)
 * and
 *
 * Copyright 2006 by Timothy E. Holy
 */

void gradEwork(const float *g[],int n_dims,const float *psi1,const int *sz,const float *psi2,const int *sz_psi2,const float *psig,const float *sqrtdetJ,float *grad[]);

#define MAX_DIMS 3

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *g[MAX_DIMS],*psi1,*psi2,*psig,*sqrtdetJ;
  float *grad[MAX_DIMS];
  int n_dims,dimIterator,dimIterator2,n_dims_tmp;
  const int *sz_tmp;
  int sz[MAX_DIMS],sz_work[MAX_DIMS],sz_psi2[MAX_DIMS];
  const mxArray *curarg,*cellarg;
  mxArray *outarg;

  if (nrhs < 4 || nrhs > 5)
    mexErrMsgTxt("register_gradE: requires four or five inputs");
  if (nlhs != 1)
    mexErrMsgTxt("register_gradE: requires one output");

  // Parse the inputs
  // g: (this is the hard one!)
  curarg = prhs[0];
  if (!mxIsCell(curarg))
    mexErrMsgTxt("register_gradE: g must be a cell array");
  n_dims = mxGetNumberOfElements(curarg);
  if (n_dims < 1 || n_dims > 3)
    mexErrMsgTxt("register_gradE: g must represent dimensionality (have size of) 1, 2, or 3");
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
    // Get the cell element
    cellarg = mxGetCell(curarg,dimIterator);
    if (cellarg == NULL)
      mexErrMsgTxt("register_gradE: error getting cell of g{:}");
    // Check the data type & size
    if (!mxIsNumeric(cellarg) || mxIsComplex(cellarg) || !mxIsSingle(cellarg))
      mexErrMsgTxt("register_gradE: element of g{:} must be a real single-precision array");
    sz_tmp = mxGetDimensions(cellarg);
    n_dims_tmp = mxGetNumberOfDimensions(cellarg);
    if (dimIterator == 0) {
      n_dims_tmp = skip_unity_dimensions(sz_tmp,n_dims_tmp,sz,MAX_DIMS);
      if (n_dims_tmp != n_dims)
	mexErrMsgTxt("register_gradE: elements of g{:} must have dimensionality  of length(g)");
    }
    else {
      n_dims_tmp = skip_unity_dimensions(sz_tmp,n_dims_tmp,sz_work,MAX_DIMS);
      if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
	mexErrMsgTxt("register_gradE: sizes of g{:} must be identical");
    }
    // Extract the data pointer
    g[dimIterator] = (float *) mxGetData(cellarg);
  }
    
  // psi1:
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_gradE: psi1 must be a real single-precision array");
  n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				     mxGetNumberOfDimensions(curarg),
				     sz_work,MAX_DIMS);
  if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
    mexErrMsgTxt("register_gradE: psi1 does not have the same dimensions as g");
  psi1 = (float *) mxGetData(curarg);

  // psi2:
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_gradE: psi2 must be a real single-precision array");
  n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				     mxGetNumberOfDimensions(curarg),
				     sz_psi2,MAX_DIMS);
  if (n_dims_tmp != n_dims)
    mexErrMsgTxt("register_gradE: psi2 does not have the same # of dimensions as g");
  psi2 = (float *) mxGetData(curarg);

  // psig:
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_gradE: psig must be a real single-precision array");
  n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				     mxGetNumberOfDimensions(curarg),
				     sz_work,MAX_DIMS);
  if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
    mexErrMsgTxt("register_gradE: psig does not have the same dimensions as g");
  psig = (float *) mxGetData(curarg);

  // sqrtdetJ:
  if (nrhs > 4) {
    curarg = prhs[4];
    if (mxIsEmpty(curarg))
      sqrtdetJ = NULL;
    else {      
      if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
	mexErrMsgTxt("register_gradE: sqrtdetJ must be a real single-precision array");
      n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
					 mxGetNumberOfDimensions(curarg),
					 sz_work,MAX_DIMS);
      if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
	mexErrMsgTxt("register_gradE: sqrtdetJ does not have the same dimensions as g");
      sqrtdetJ = (float *) mxGetData(curarg);
    }
  } else
    sqrtdetJ = NULL;
    

  // Set up the output cell array
  plhs[0] = mxCreateCellMatrix(1,n_dims);
  if (plhs[0] == NULL)
    mexErrMsgTxt("register_gradE: error creating output (grad{:})");
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
    outarg = mxCreateNumericArray(n_dims,sz,mxSINGLE_CLASS,mxREAL);
    if (outarg == NULL)
      mexErrMsgTxt("register_gradE: error creating element of grad{:}");
    mxSetCell(plhs[0],dimIterator,outarg);
    grad[dimIterator] = (float *) mxGetData(outarg);
  }

  // Do the actual work
  gradEwork(g,n_dims,psi1,sz,psi2,sz_psi2,psig,sqrtdetJ,grad);

  return;
}


/*
 * Do the actual computation.
 * It uses explicit formulas for inverse-transpose in dimensions 1, 2, and 3. 
 * Presumably the overhead for LU-decomposition would not be justifiable.
 */
void gradEwork(const float *g[],int n_dims,const float *psi1,const int *sz,const float *psi2,const int *sz_psi2,const float *psig,const float *sqrtdetJ,float *grad[])
{
  bool covariant;
/*  int *pixel_skip,*bit_value;*/
  int pixel_skip[MAX_DIMS+1],bit_value[MAX_DIMS+1];
  int dimIterator,dimIterator2,interpIterator,n_interp_terms,pixel_offset,gxf,tmpskip;
/*  float *J, *Jti, *alpha;*/
  float J[MAX_DIMS*MAX_DIMS], Jti[MAX_DIMS*MAX_DIMS], alpha[MAX_DIMS];
  const float *psi1_end, *psi2_tmp;
  float coef,tmp,interp_sum,interp_coef,coef_alt,detJ;

  covariant = (sqrtdetJ != NULL);
/*
//   J = (float *) mxMalloc(n_dims*n_dims*sizeof(float));
//   Jti = (float *) mxMalloc(n_dims*n_dims*sizeof(float));
//   alpha = (float *) mxMalloc(n_dims*sizeof(float));
//   pixel_skip = (int *) mxMalloc((n_dims+1)*sizeof(int));
//   bit_value = (int *) mxMalloc((n_dims+1)*sizeof(int));
*/

  pixel_skip[0] = 1;
  bit_value[0] = 1;
  for (dimIterator = 1; dimIterator <= n_dims; dimIterator++) {
    pixel_skip[dimIterator] = pixel_skip[dimIterator-1]*sz[dimIterator-1];
    bit_value[dimIterator] = 1 << dimIterator;
  }
  psi1_end = psi1+pixel_skip[n_dims];
  n_interp_terms = bit_value[n_dims];
  //
  // Loop over output pixels (would start multithreading here)
  // (Note that the incrementing of g[] and grad[] is occuring inside
  // the loop rather than in the conventional spot in the "for"
  // construct.)
  //
  pixIterator pixI(sz,n_dims);
  for (; psi1 < psi1_end; psi1++, psig++, sqrtdetJ++, pixI++) {
  //for (; psi1 < psi1_end; psi1++, psig++, sqrtdetJ++) {
    coef = 2*(*psig - *psi1);
    // Check to see whether warped image data is even available for this pixel (e.g., if the warp makes part of the image go "offscreen"); this also checks whether we're on the edge of the image. This also obviates the need to worry about whether nearest-neighbors are available.
    if (mxIsNaN(coef)) {
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
	*(grad[dimIterator]) = coef;
        if (covariant) {
          // Have to propagate NaNs inward a bit to compensate for edge problems
          tmpskip = pixel_skip[dimIterator];
          if (pixI.coord(dimIterator) > 0)
	   grad[dimIterator][-tmpskip] = coef;
          if (pixI.coord(dimIterator) < sz[dimIterator]-1)
	   grad[dimIterator][tmpskip] = coef;
        }
        grad[dimIterator]++;
	g[dimIterator]++;
      }
      continue;  // go on to the next pixel
    }
    if (covariant) {
      // Calculate the derivative of the detJ term.  This requires us
      // to recalculate J; however, the overhead of this is very small
      // compared to the total calculation, and if we were to instead
      // store J from computing psig, it would add tremendously to the
      // memory requirements.
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
	for (dimIterator2 = 0; dimIterator2 < n_dims; dimIterator2++)
	  J[dimIterator+dimIterator2*n_dims] = 
	    (g[dimIterator][pixel_skip[dimIterator2]] 
	     - g[dimIterator][-pixel_skip[dimIterator2]])/2;
      // Now compute the components of inv(J).  We do this by explicit
      // formulas for dimensions 1,2,and 3, this presumably being
      // faster than a general algorithm like LU decomposition. Note
      // that we recompute the determinant, rather than squaring
      // sqrtdetJ, because of the absolute value: sqrtdetJ is really
      // sqrt(abs(det(J))).
      if (n_dims == 1) {
	Jti[0] = 1/J[0];
	//tmp = coef * *psig * Jti[0]/4;
	//grad[0][1] += tmp;
	//grad[0][-1] -= tmp;
      }
      else if (n_dims == 2) {
	detJ = J[0]*J[3]-J[1]*J[2];
	Jti[0] = J[3] / detJ;
	Jti[1] = -J[2] / detJ;
	Jti[2] = -J[1] / detJ;
	Jti[3] = J[0] / detJ;
      }
      else {
	Jti[0] = (J[4]*J[8]-J[7]*J[5]);
	Jti[3] = -(J[1]*J[8]-J[7]*J[2]);
	Jti[6] = (J[1]*J[5]-J[4]*J[2]);
	detJ = J[0]*Jti[0] + J[3]*Jti[3] + J[6]*Jti[6];
	Jti[0] /= detJ;
	Jti[3] /= detJ;
	Jti[6] /= detJ;
	Jti[1] = -(J[3]*J[8]-J[6]*J[5]) / detJ;
	Jti[4] = (J[0]*J[8]-J[6]*J[2]) / detJ;
	Jti[7] = -(J[0]*J[5]-J[3]*J[2]) / detJ;
	Jti[2] = (J[3]*J[7]-J[6]*J[4]) / detJ;
	Jti[5] = -(J[0]*J[7]-J[6]*J[1]) / detJ;
	Jti[8] = (J[0]*J[4]-J[3]*J[1]) / detJ;
      }
      coef_alt = coef * *psig / 4;
      //coef_alt = *sqrtdetJ/4;  // For debugging: energy = sum sqrt(detJ)
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
	for (dimIterator2 = 0; dimIterator2 < n_dims; dimIterator2++) {
	  tmp = coef_alt * Jti[dimIterator + dimIterator2*n_dims];
          tmpskip = pixel_skip[dimIterator2];
          grad[dimIterator][tmpskip] += tmp;
          grad[dimIterator][-tmpskip] -= tmp;
	}
      // Now we're done with the J terms. One last step is we can
      // adjust the coefficient in perparation for the interpolation
      // component (since this, too, is only done for the covariant
      // theory)
      coef *= *sqrtdetJ;
    }  // if(covariant)

    // Extract the floor and fractional components of g(x), and set
    // psi2_tmp to point to psi2(floor(g(x))).  This is necessary
    // preparation for the interpolation component.
    psi2_tmp = psi2;
    for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
      gxf = (int) floorf(*(g[dimIterator]));
      psi2_tmp += (gxf-1) * pixel_skip[dimIterator];   // The -1 is because MATLAB is unit-offset
      alpha[dimIterator] = *(g[dimIterator]) - gxf;
    }
    // Loop over gradient dimensions
    for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
      // Compute the gradient for the interpolation term. This is the
      // single slowest step (requiring the most flops) because there
      // are d+1 multiplies in each of 2^d terms
      interp_sum = 0;
      for (interpIterator = 0; interpIterator < n_interp_terms; interpIterator++) {
	interp_coef = 1;
	pixel_offset = 0;
	for (dimIterator2 = 0; dimIterator2 < n_dims; dimIterator2++) {
	  if (dimIterator2 == dimIterator)
	    interp_coef *= (2*(interpIterator&bit_value[dimIterator2])-1);
	  else {
	    // Which is faster, an if or the multiply?
	    if (interpIterator&bit_value[dimIterator2])
	      interp_coef *= alpha[dimIterator2];
	    else
	      interp_coef *= 1-alpha[dimIterator2];
	  }
	  if (interpIterator&bit_value[dimIterator2])
	    pixel_offset += pixel_skip[dimIterator2];
	}
	interp_sum += interp_coef * psi2_tmp[pixel_offset];
      }
      *(grad[dimIterator]) += coef*interp_sum;
    }
    
    // Increment the grad[] and g[] pointers
    for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
      grad[dimIterator]++;
      g[dimIterator]++;
    }
  } // loop over pixels
/*  mxFree(J);
  mxFree(Jti);
  mxFree(alpha);
  mxFree(pixel_skip);
  mxFree(bit_value);*/
}

#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 5;
  const int n_outputs = 1;
  const mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "g",
    "psi1",
    "psi2",
    "psig",
    "sqrtdetJ"
  };
  const char *output_names[] = {
    "grad"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables g, psi1, psi2, psig, and sqrtdetJ\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  printf("Output file just before calling mexfcn: %s\n",fileout);
  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,input);

  // Save the outputs
  printf("Output file just before save: %s\n",fileout);
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
