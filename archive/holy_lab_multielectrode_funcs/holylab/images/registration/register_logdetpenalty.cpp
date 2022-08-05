#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#else
#include "mex.h"
#endif

#include <math.h>   /* for log */
#include "imiterators.cxx"  /* for pixIterator */
#include "image_utilities.h"  /* for skip_unity_dimensions */

/*
 * register_logdetpenalty: compute the gradient of the image difference error
 * with respect to coordinate deformation
 *
 * Syntax:
 *   [detJnew,val,grad,hessi,hessj,hessv] = register_logdetpenalty(u,c,detJprev)
 *   [detJnew,val,grad,hessi,hessv] = register_logdetpenalty(u,c,detJprev)
 *      (this second syntax computes only the diagonal elements of the hessian)
 * where
 *   detJprev is optional, c defaults to 1
 * and
 *
 * Copyright 2010 by Timothy E. Holy
 */


#define MAX_DIMS 3
// The next should be 2^MAX_DIMS
#define MAX_NBRS 8


// Function declarations
template <class dataType>
void mexWrapper(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[],
		 mxClassID classID,dataType (*logfunc)(dataType));
template <class dataType>
void reg_ldp_work(int n_dims,const int szu[],const dataType *u,dataType c,const dataType *detJprev,dataType *detJnew,dataType &val,dataType *grad,dataType *hessi,dataType *hessj,dataType *hessv,dataType (*logfunc)(dataType),dataType infval);



/*
 * This is the Matlab "outer" wrapper, that initiates templates depending on datatype
 */
void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  mxClassID classID;

  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt("register_logdetpenalty: requires one to three inputs");
  if ((nlhs > 3 && nlhs < 5) || nlhs > 6)
    mexErrMsgTxt("register_logdetpenalty: must have 0-3, 5, or 6 outputs");
  if (nlhs == 0)
    return;

  // Determine the data type
  curarg = prhs[0];
  if (!mxIsNumeric(curarg))
    mexErrMsgTxt("register_logdetpenalty: u must be numeric");
  classID = mxGetClassID(curarg);

  // Check data type consistency
  if (nrhs > 2) {
    curarg = prhs[2];
    if (!mxIsNumeric(curarg) || mxGetClassID(curarg) != classID)
      mexErrMsgTxt("register_logdetpenalty: the data types of u and detJprev must agree");
  }

  if (classID == mxSINGLE_CLASS)
    mexWrapper<float>(nlhs,plhs,nrhs,prhs,mxSINGLE_CLASS,logf);
  else if (classID == mxDOUBLE_CLASS)
    mexWrapper<double>(nlhs,plhs,nrhs,prhs,mxDOUBLE_CLASS,log);
  else
    mexErrMsgTxt("register_logdetpenalty: supports only single and double");
}


  /*
   * This is the Matlab "inner" wrapper
   */
template <class dataType>
void mexWrapper(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[],
		 mxClassID classID,dataType (*logfunc)(dataType))
{
  int n_dims,dimIndex,n_dims_tmp;
  const dataType *u,*detJprev;
  dataType *detJnew,*grad,*hessi,*hessj,*hessv,*valP;
  dataType c,val;
  int sz[MAX_DIMS+1];
  int szu[MAX_DIMS+1];
  int dimList[MAX_DIMS+1];
  const mxArray *curarg;
  const mwSize *sztmp;

  // Parse the arguments
  // u
  curarg = prhs[0];
  n_dims = mxGetNumberOfDimensions(curarg);
  if (n_dims > MAX_DIMS+1)
    mexErrMsgTxt("register_logdetpenalty: u has more than the maximum number of dimensions (fix with recompile?)");
  sztmp = mxGetDimensions(curarg);
  // By getting rid of unity dimensions now, we ensure that in the
  // work function that the pixIterator for detJ has the same
  // dimensionality that the pixIterator for u has (we would otherwise
  // have trouble with dimensions of size 2).
  n_dims = skip_unity_dimensions_index(sztmp,n_dims,dimList,MAX_DIMS+1)-1;
  if (n_dims <= 0) {
    n_dims = 1;
    dimList[1] = 1;
  }
  for (dimIndex = 0; dimIndex <= n_dims; dimIndex++)
    szu[dimIndex] = sztmp[dimList[dimIndex]];
  if (szu[n_dims] != n_dims)
    mexErrMsgTxt("register_logdetpenalty: the shape of u is inconsistent with the dimensionality");
  if (n_dims < 1 || n_dims > MAX_DIMS)
    mexErrMsgTxt("register_logdetpenalty: u must represent dimensionality 1, 2, or 3");
  u = (dataType*) mxGetData(curarg);

  // c
  if (nrhs > 1)
    c = mxGetScalar(prhs[1]);
  else
    c = 1;

  // detJprev
  detJprev = NULL;
  if (nrhs > 2) {
    curarg = prhs[2];
    //if (mxGetNumberOfDimensions(curarg) != n_dims)
    //  mexErrMsgTxt("register_logdetpenalty: detJprev does not have the right number of dimensions");
    n_dims_tmp = mxGetNumberOfDimensions(curarg);
    if (n_dims_tmp > n_dims)
      mexErrMsgTxt("register_logdetpenalty: detJprev has too many dimensions");
    if (n_dims_tmp < n_dims) {
      // If the last dimensions are all 1, then it's actually OK
      for (dimIndex = n_dims_tmp; dimIndex < n_dims; dimIndex++)
	if (szu[dimIndex] > 2)
	  mexErrMsgTxt("register_logdetpenalty: detJprev has too few dimensions");
    }
    sztmp = mxGetDimensions(curarg);
    for (dimIndex = 0; dimIndex < n_dims_tmp; dimIndex++)
      if (sztmp[dimList[dimIndex]] != szu[dimIndex]-1)
	mexErrMsgTxt("register_logdetpenalty: the size of detJprev is not consistent with u");
    detJprev = (dataType*) mxGetData(curarg);
  }

  // Allocate memory for outputs
  detJnew = NULL;
  valP = NULL;
  grad = NULL;
  hessi = hessj = hessv = NULL;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    sz[dimIndex] = szu[dimIndex]-1;
  // detJnew
  if (nlhs > 0) {
    plhs[0] = mxCreateNumericArray(n_dims,sz,classID,mxREAL);
    detJnew = (dataType *) mxGetData(plhs[0]);
  }
  // val
  if (nlhs > 1) {
    plhs[1] = mxCreateNumericMatrix(1,1,classID,mxREAL);
    valP = (dataType *) mxGetData(plhs[1]);
  }
  // grad
  if (nlhs > 2) {
    plhs[2] = mxCreateNumericArray(n_dims+1,szu,classID,mxREAL);
    grad = (dataType *) mxGetData(plhs[2]);
  }
  // hess
  if (nlhs > 3) {
    int n_nbrs = 1;
    long n_elem_detJ = 1;
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
      n_nbrs *= 2;
      n_elem_detJ *= sz[dimIndex];
    }
    if (nlhs == 6) {
      long n_entries_hess = n_nbrs*n_nbrs*n_elem_detJ*n_dims*n_dims;
      //mexPrintf("n_entries_hess = %d, n_elem_detJ %d, n_nbrs %d\n",n_entries_hess,n_elem_detJ,n_nbrs);
      plhs[3] = mxCreateNumericMatrix(n_entries_hess,1,classID,mxREAL);
      hessi = (dataType *) mxGetData(plhs[3]);
      plhs[4] = mxCreateNumericMatrix(n_entries_hess,1,classID,mxREAL);
      hessj = (dataType *) mxGetData(plhs[4]);
      plhs[5] = mxCreateNumericMatrix(n_entries_hess,1,classID,mxREAL);
      hessv = (dataType *) mxGetData(plhs[5]);
    } else {
      long n_entries_hess = n_nbrs*n_elem_detJ*n_dims;
      //mexPrintf("n_entries_hess = %d, n_elem_detJ %d, n_nbrs %d\n",n_entries_hess,n_elem_detJ,n_nbrs);
      plhs[3] = mxCreateNumericMatrix(n_entries_hess,1,classID,mxREAL);
      hessi = (dataType *) mxGetData(plhs[3]);
      plhs[4] = mxCreateNumericMatrix(n_entries_hess,1,classID,mxREAL);
      hessv = (dataType *) mxGetData(plhs[4]);
    }
  }

  // Execute the function
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    sz[dimIndex] = szu[dimIndex];
  reg_ldp_work(n_dims,sz,u,c,detJprev,detJnew,val,grad,hessi,hessj,hessv,logfunc,(dataType) mxGetInf());
  if (nlhs > 1)
    *valP = val;
}

/*
 * Do the actual computation.
 * It uses explicit formulas for inverse-transpose in dimensions 1, 2, and 3. 
 * Presumably the overhead for LU-decomposition would not be justifiable.
 */
template <class dataType>
void reg_ldp_work(int n_dims,const int szu[],const dataType *u,dataType c,const dataType *detJprev,dataType *detJnew,dataType &val,dataType *grad,dataType *hessi,dataType *hessj,dataType *hessv,dataType (*logfunc)(dataType),dataType infval)
{
  size_t umemoffset[MAX_NBRS];
  int coordoffset[MAX_NBRS][MAX_DIMS];
  int sztmp[MAX_DIMS];
  size_t uIndex,skip,skip2;
  int n_nbrs,n_nbrs_half;
  int tmpI,dimIndex,dimIndex1,dimIndex2,dimIndex3,nbrIndex,nbrIndex2,Jindex,Jindex2;
  dataType J[MAX_DIMS*MAX_DIMS], Jti[MAX_DIMS*MAX_DIMS];
  dataType L,detJ;
  dataType coef2,coef4;
  dataType tmp,tmp1,tmp2;
  dataType *hessiP;

  hessiP = hessi;
  pixIterator pIu(szu,n_dims,false);
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    sztmp[dimIndex] = szu[dimIndex]-1;
  pixIterator pIdet(sztmp,n_dims,false);// see comment about pixIterator above
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    sztmp[dimIndex] = 2;
  pixIterator offset(sztmp,n_dims);
  
  // Pre-calculate the various memory offsets from neighbors
  n_nbrs = 1;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    n_nbrs *= 2;
  n_nbrs_half = n_nbrs/2;
  for (nbrIndex = 0, offset.restart(); nbrIndex < n_nbrs; nbrIndex++,offset++) {
    tmpI = 0;
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
      tmpI += offset.coord(dimIndex)*pIu.dimSkip(dimIndex);
      coordoffset[nbrIndex][dimIndex] = offset.coord(dimIndex);
    }
    umemoffset[nbrIndex] = tmpI;
  }
  coef2 = 2.0/n_nbrs;  // 1/2^(d-1)
  coef4 = 2*coef2;

  // Initialize
  val = 0;

  // Iterate over pixels of determinant
  for (; !pIdet.at_end(); pIdet++) {
    // Get the index for u that corresponds to the current half-grid point
    uIndex = pIu.offset_from_coords(pIdet);
    // Calculate the Jacobian at the current point
    for (skip = 0, dimIndex1 = 0; dimIndex1 < n_dims; dimIndex1++, skip+=pIu.numel())
      for (dimIndex2 = 0,Jindex=dimIndex1; dimIndex2 < n_dims; dimIndex2++,Jindex+=n_dims) {
	tmp = 0;
	for (nbrIndex = 0; nbrIndex < n_nbrs; nbrIndex++) {
	  tmp1 = u[skip+uIndex+umemoffset[nbrIndex]];
	  if (coordoffset[nbrIndex][dimIndex2])
	    tmp += tmp1;
	  else
	    tmp -= tmp1;
	}
	if (dimIndex1 == dimIndex2)
	  J[Jindex] = tmp*coef2 + 1;  // add identity for u->g
	else
	  J[Jindex] = tmp*coef2;
      }

    // Compute detJ
    if (n_dims == 1) {
      detJ = J[0];
    }
    else if (n_dims == 2) {
      detJ = J[0]*J[3]-J[1]*J[2];
    }
    else {
      // We'll get a start on Jti, where appropriate, as we compute detJ
      Jti[0] = (J[4]*J[8]-J[7]*J[5]);
      Jti[3] = -(J[1]*J[8]-J[7]*J[2]);
      Jti[6] = (J[1]*J[5]-J[4]*J[2]);
      detJ = J[0]*Jti[0] + J[3]*Jti[3] + J[6]*Jti[6];
    }
    detJnew[pIdet] = detJ;

    // Compute contribution to penalty value
    L = detJ/c;
    if (detJprev != NULL)
      L *= detJprev[pIdet];
    if (L < 0)
      L = infval;
    else
      L = logfunc(L);
    val += L*L;

    // Compute gradient
    if (grad != NULL) {
      // Compute Jti (J^T^(-1))
      if (n_dims == 1) {
	Jti[0] = 1.0/detJ;
      }
      else if (n_dims == 2) {
	Jti[0] = J[3] / detJ;
	Jti[1] = -J[2] / detJ;
	Jti[2] = -J[1] / detJ;
	Jti[3] = J[0] / detJ;
      }
      else {
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
      // Compute contributions to grad
      L *= coef4;
      for (nbrIndex = 0; nbrIndex < n_nbrs_half; nbrIndex++) {
	for (dimIndex = 0, skip=0; dimIndex < n_dims; dimIndex++,skip+=pIu.numel()) {
	  tmp = 0;
	  for (dimIndex2 = 0,Jindex=dimIndex; dimIndex2 < n_dims; dimIndex2++,Jindex+=n_dims) {
	    if (coordoffset[nbrIndex][dimIndex2])
	      tmp += Jti[Jindex];
	    else
	      tmp -= Jti[Jindex];
	  }
	  tmp *= L;
	  grad[uIndex+skip+umemoffset[nbrIndex]] += tmp;
	  grad[uIndex+skip+umemoffset[n_nbrs-nbrIndex-1]] -= tmp;
	}
      }
      if (hessi != NULL) {
	dataType coefhess = 2.0*coef2*coef2/pIdet.numel();
	for (nbrIndex = 0; nbrIndex < n_nbrs_half; nbrIndex++)
	  for (nbrIndex2 = 0; nbrIndex2 < n_nbrs_half; nbrIndex2++) {
	    if (hessj == NULL && nbrIndex != nbrIndex2)
	      continue;
	    for (dimIndex = 0, skip = 0; dimIndex < n_dims; dimIndex++, skip+=pIu.numel())
	      for (dimIndex2 = 0, skip2 = 0; dimIndex2 < n_dims; dimIndex2++, skip2+=pIu.numel()) {
		if (hessj == NULL && dimIndex != dimIndex2)
		  continue;
		tmp = tmp2 = 0;
		for (dimIndex3 = 0, Jindex=dimIndex,Jindex2=dimIndex2; dimIndex3 < n_dims; dimIndex3++,Jindex+=n_dims,Jindex2+=n_dims) {
		  if (coordoffset[nbrIndex][dimIndex3])
		    tmp += Jti[Jindex];
		  else
		    tmp -= Jti[Jindex];
		  if (coordoffset[nbrIndex2][dimIndex3])
		    tmp2 += Jti[Jindex2];
		  else
		    tmp2 -= Jti[Jindex2];
		}
		dataType prod = tmp*tmp2*coefhess;
		long c1 = uIndex+skip+umemoffset[nbrIndex];
		if (hessj != NULL) {
		  long c2 = uIndex+skip2+umemoffset[nbrIndex2];
		  *hessi++ = c1; *hessj++ = c2; *hessv++ = prod;
		  c2 = uIndex+skip2+umemoffset[n_nbrs-nbrIndex2-1];
		  *hessi++ = c1; *hessj++ = c2; *hessv++ = -prod;
		  c1 = uIndex+skip+umemoffset[n_nbrs-nbrIndex-1];
		  *hessi++ = c1; *hessj++ = c2; *hessv++ = prod;
		  c2 = uIndex+skip2+umemoffset[nbrIndex2];
		  *hessi++ = c1; *hessj++ = c2; *hessv++ = -prod;
		} else {
		  *hessi++ = c1; *hessv++ = prod;
		  c1 = uIndex+skip+umemoffset[n_nbrs-nbrIndex-1];
		  *hessi++ = c1; *hessv++ = prod;
		}
		//mexPrintf("hessi: %d\n",hessi-hessiP);
	      }
	  }
      }
    }
  }
  // Divide by the number of voxels in the determinant (to normalize)
  val /= pIdet.numel();
  if (grad != NULL) {
    for (uIndex = 0; uIndex < pIu.numel()*n_dims; uIndex++)
      grad[uIndex] /= pIdet.numel();
  }
  // no need to do hess, it was fixed above
}

#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 3;
  const int n_outputs = 6;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "u",
    "c",
    "detJprev"
  };
  const char *output_names[] = {
    "detJnew",
    "val",
    "grad",
    "hessi",
    "hessj",
    "hessv"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables u, c, and detJprev\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray **) input);

  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif

