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
 * register_gradE: compute the gradient of the registration error with
 * respect to coordinate deformation. This includes the regularization
 * components.
 *
 * Syntax:
 *   grad = register_gradE(g,psi1,psi2,psig,weight,deltax,lambda,sqrtdetJ,covtflag)
 * where
 *  weight is weight for each pixel in psi1/psig
 *  deltax is spacing along each coordinate
 *  lambda = [lambda1 lambda2] (lambda1 sets penalty on detJ being different from 1, lambda2 sets penalty on grad(sqrtdetJ))
 * and
 *
 * Copyright 2006-2007 by Timothy E. Holy
 */

void gradEwork(const float *g[],int n_dims,const float *psi1,const int *sz,const float *psi2,const int *sz_psi2,const float *psig,const float *weight,const double *deltax,double lambda1,double lambda2,const float *sqrtdetJ,bool covariant,float *grad[]);

#define MAX_DIMS 3

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *g[MAX_DIMS],*psi1,*psi2,*psig,*sqrtdetJ, *weight;
  const double *deltax, *lambda_ptr;
  float *grad[MAX_DIMS];
  double lambda1,lambda2;
  int n_dims,dimIterator,dimIterator2,n_dims_tmp;
  const int *sz_tmp;
  int sz[MAX_DIMS],sz_work[MAX_DIMS],sz_psi2[MAX_DIMS];
  bool covariant;
  const mxArray *curarg,*cellarg;
  mxArray *outarg;

  if (nrhs != 9)
    mexErrMsgTxt("register_gradE: requires nine inputs");
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

  
  // weight:
  curarg = prhs[4];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_gradE: mask must be a real single-precision array");
  n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				     mxGetNumberOfDimensions(curarg),
				     sz_work,MAX_DIMS);
  if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
    mexErrMsgTxt("register_gradE: mask does not have the same dimensions as g");
  weight = (float *) mxGetData(curarg);


  // deltax
  curarg = prhs[5];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != n_dims)
    mexErrMsgTxt("register_gradE: deltax must have length equal to the length of g");
  deltax = mxGetPr(curarg);
  

  // lambda
  curarg = prhs[6];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg) || mxGetNumberOfElements(curarg) != 2)
    mexErrMsgTxt("register_gradE: lambda must be a real double-precision 2-vector");
  lambda_ptr = mxGetPr(curarg);
  lambda1 = lambda_ptr[0];
  lambda2 = lambda_ptr[1];
  

  // sqrtdetJ:
  curarg = prhs[7];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("register_gradE: sqrtdetJ must be a real single-precision array");
  n_dims_tmp = skip_unity_dimensions(mxGetDimensions(curarg),
				     mxGetNumberOfDimensions(curarg),
				     sz_work,MAX_DIMS);
  if (!validate_dimensions(sz,n_dims,sz_work,n_dims_tmp))
    mexErrMsgTxt("register_gradE: sqrtdetJ does not have the same dimensions as g");
  sqrtdetJ = (float *) mxGetData(curarg);
    

  // covtflag
  curarg = prhs[8];
  covariant = (bool) mxGetScalar(curarg);


  // Set up the output cell array
  plhs[0] = mxCreateCellMatrix(1,n_dims);
  if (plhs[0] == NULL)
    mexErrMsgTxt("register_gradE: error creating output (grad{:})");
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
    outarg = mxCreateNumericArray(mxGetNumberOfDimensions(cellarg),
				  mxGetDimensions(cellarg),
				  mxSINGLE_CLASS,mxREAL);
    if (outarg == NULL)
      mexErrMsgTxt("register_gradE: error creating element of grad{:}");
    mxSetCell(plhs[0],dimIterator,outarg);
    grad[dimIterator] = (float *) mxGetData(outarg);
  }

  // Do the actual work
  gradEwork(g,n_dims,psi1,sz,psi2,sz_psi2,psig,weight,deltax,lambda1,lambda2,sqrtdetJ,covariant,grad);

  return;
}


/*
 * Do the actual computation.  It uses explicit formulas for
 * inverse-transpose in dimensions 1, 2, and 3.  Presumably the
 * overhead for LU-decomposition would not be justifiable; hence the
 * restriction to these dimensions.
 */
void gradEwork(const float *g[],int n_dims,const float *psi1,const int *sz,const float *psi2,const int *sz_psi2,const float *psig,const float *weight,const double *deltax,double lambda1,double lambda2,const float *sqrtdetJ,bool covariant,float *grad[])
{
/*  int *pixel_skip,*bit_value;*/
  int pixel_skip[MAX_DIMS+1],bit_value[MAX_DIMS+1];
  int dimIterator,dimIterator2,interpIterator,n_interp_terms,pixel_offset,gxf,tmpskip;
/*  float *J, *Jti, *alpha;*/
  float J[MAX_DIMS*MAX_DIMS], Jti[MAX_DIMS*MAX_DIMS], alpha[MAX_DIMS];
  const float *psi1_end, *psi2_tmp;
  float img_diff,tmp,interp_sum,interp_coef,coef,detJ,xi_lap,denom;
  bool in_interior;

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
  for (; psi1 < psi1_end; psi1++, psig++, weight++, sqrtdetJ++, pixI++) {
    if (*weight == 0)
      img_diff = 0;   // set to 0 even if psig is NaN
    else
      img_diff = (*psig - *psi1) * *weight;
    if (*weight > 0) {
      //
      // The derivative with respect to g of the interpolated pixel values
      //
      // Extract the floor and fractional components of g(x), and set
      // psi2_tmp to point to psi2(floor(g(x))).  This is necessary
      // preparation for the interpolation component.
      psi2_tmp = psi2;
      in_interior = true;
      for (dimIterator = 0; dimIterator < n_dims && in_interior; dimIterator++) {
	gxf = (int) floorf(*(g[dimIterator]));
	// The following lines are influenced by Matlab's unit-offset
	if (gxf < 1 || gxf >= sz_psi2[dimIterator])
	  in_interior = false;
	psi2_tmp += (gxf-1) * pixel_skip[dimIterator];
	alpha[dimIterator] = *(g[dimIterator]) - gxf;
      }
      if (in_interior) {
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
	  if (covariant)
	    interp_sum *= *sqrtdetJ;
	  *(grad[dimIterator]) += 2*img_diff*interp_sum;
	}
      }
    }
    if ((covariant || lambda1 > 0 || lambda2 > 0) && !pixI.on_edge()) {
      // Calculate the derivative of the detJ term.  This requires us
      // to recalculate J (it was calculated for psig); however, the
      // overhead of this is very small compared to the total
      // calculation, and if we were to instead store J from computing
      // psig, it would add tremendously to the memory requirements.
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
	for (dimIterator2 = 0; dimIterator2 < n_dims; dimIterator2++) {
          tmpskip = pixel_skip[dimIterator2];
	  if (pixI.coord(dimIterator2) > 0 &&
	      pixI.coord(dimIterator2) < sz[dimIterator2]-1)
	    J[dimIterator+dimIterator2*n_dims] = 
	      (g[dimIterator][tmpskip] - g[dimIterator][-tmpskip])/2;
	  else if (pixI.coord(dimIterator2) == 0)
	    J[dimIterator+dimIterator2*n_dims] = 
	      g[dimIterator][tmpskip] - *(g[dimIterator]);
	  else
	    J[dimIterator+dimIterator2*n_dims] = 
	      *(g[dimIterator]) - g[dimIterator][-tmpskip];
	}
	    
      // Now compute the components of inv(J).  We do this by explicit
      // formulas for dimensions 1,2,and 3, this presumably being
      // faster than a general algorithm like LU decomposition. Note
      // that we recompute the determinant, rather than squaring
      // sqrtdetJ, because of the absolute value: sqrtdetJ is really
      // sqrt(abs(det(J))).
      // Note furthermore: because of composition, it might not be
      // true that sqrtdetJ = sqrt(abs(detJ))! So we definitely need
      // to calculate the determinant again.
      if (n_dims == 1) {
	Jti[0] = 1/J[0];
	//tmp = img_diff * *psig * Jti[0]/4;
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
      // Calculate the laplacian of sqrtdetJ = xi
      xi_lap = 0;
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
	denom = deltax[dimIterator];
	denom = denom*denom;
	if (pixI.coord(dimIterator) > 0)
	  xi_lap += (sqrtdetJ[-pixel_skip[dimIterator]] - *sqrtdetJ)/denom;
	if (pixI.coord(dimIterator) < sz[dimIterator]-1)
	  xi_lap += (sqrtdetJ[pixel_skip[dimIterator]] - *sqrtdetJ)/denom;
      }
      coef = *sqrtdetJ * (lambda1 * (*sqrtdetJ-1.0) - lambda2 * xi_lap)/2;
      if (*weight > 0)    // to prevent NaNs
	coef += img_diff * *psig / 2;
      //coef = *sqrtdetJ/4;  // For debugging: energy = sum sqrt(detJ)
      for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
	for (dimIterator2 = 0; dimIterator2 < n_dims; dimIterator2++) {
	  tmp = coef * Jti[dimIterator + dimIterator2*n_dims]/deltax[dimIterator2];
          tmpskip = pixel_skip[dimIterator2];
	  if (pixI.coord(dimIterator2) > 0 && 
	      pixI.coord(dimIterator2) < sz[dimIterator2]-1) {
	    grad[dimIterator][tmpskip] += tmp;
	    grad[dimIterator][-tmpskip] -= tmp;
	  } else if (pixI.coord(dimIterator2) == 0) {
	    grad[dimIterator][tmpskip] += 2*tmp;
	    *(grad[dimIterator]) -= 2*tmp;
	  }
	  else {
	    *(grad[dimIterator]) += 2*tmp;
	    grad[dimIterator][-tmpskip] -= 2*tmp;
	  }
	}
    }  // if(covariant)

    
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
  const int n_inputs = 9;
  const int n_outputs = 1;
  const mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "g",
    "psi1",
    "psi2",
    "psig",
    "weight",
    "deltax",
    "lambda",
    "sqrtdetJ",
    "covariantflag"
  };
  const char *output_names[] = {
    "grad"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables ",argv[0]);
    for (int i = 0; i < n_inputs; i++)
      printf("%s, ",input_names[i]);
    printf("\n");
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,input);

  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
