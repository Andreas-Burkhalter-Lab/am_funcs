/*
 * See register_logdetpenalty_composition.m for syntax.
 *
 * Copyright 2012 by Timothy E. Holy and Jian Wang
 */

#ifdef MAIN 
/* Header for debugging purpose, independent from MATLAB program */
#include <stdio.h>
#include "mat.h"
#include "mat_variables.h"

#else
/* Header for using within MATLAB program */
#include "mex.h"

#endif

#include "matrix.h"
#include <stdlib.h>
#include <cmath> 						/* for log */
#include <limits> 					/* for infinity */
#include "imiterators.cxx" 	/* for pixIterator */
//#include "image_utilities.h" /* for skip_unity_dimensions */
#include "jacobianu.cpp" 		/* for Jacobian matrix related computation */
#include "qinterp.cpp" 			/* for quadratic interpolation */

#define MAX_DIMS 3
#define MAX_NBRS 8 /* equal to 2^MAX_DIMS */


/*
 * Function pre-declarations
 */

/* Inner-wrapper function called by the outer wrapper mexFunction to
   parse arguments */
template <class T>
void mexWrapper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],
		mxClassID classID);

/* Real function to do the actual computation */
template <class T, int n_dims>
void regldpc_work(const int szu[], const T *uold[],	const T *unew[], T c,
		T *ucomp[], T *LDetJ, T *valp, T *grad[], T *ucompgrad[]);

/* Check the validity of input arguments */
void check_cell_args(const mxArray *curarg, int d, mxClassID classID,
		const char *vname) {
	
	int dimIndex;
	const mxArray *newarg;

	if (!mxIsCell(curarg)) {
		mexErrMsgIdAndTxt("register:penalty", 
				"register_logdetpenalty_composition: %s must be a cell array", vname);
	}
	if (mxGetNumberOfElements(curarg) != d) {
		mexErrMsgIdAndTxt("register:penalty",	
				"register_logdetpenalty_composition: %s must have %d elements", 
				vname, d);
	}
	for (dimIndex = 0; dimIndex < d; ++dimIndex) {
		newarg = mxGetCell(curarg, dimIndex);
		if (mxGetClassID(newarg) != classID) {
			mexErrMsgIdAndTxt("register:penalty", 
					"register_logdetpenalty_composition: the elements of %s must "
					"all be of the same type as unewc{1}.",	vname);
		}
	}
}

/*
 * This is the MATLAB "outer" wrapper that initiates templates
 * depending on T
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	const mxArray *curarg;
	mxClassID classID;
	int d;

	/* Check the input & output argument numbers */
	if (nrhs < 2 || nrhs > 3) {
		mexErrMsgTxt(
				"register_logdetpenalty_composition: requires two to three inputs");
	}
	if (nlhs < 1 || nlhs == 2 || nlhs > 4) {
		mexErrMsgTxt(
				"register_logdetpenalty_composition: must have 1, 3, or 4 outputs");
	}
	if (nlhs == 0) {
		return;
	}

	/* Determine the data type of prhs[1] */
	curarg = prhs[1];
	if (!mxIsCell(curarg)) {
		mexErrMsgTxt(
				"register_logdetpenalty_composition: unewc must be a cell array");
	}
	d = mxGetNumberOfElements(curarg); /* number of cells in prhs[1] */
	curarg = mxGetCell(curarg, 0);
	if (!mxIsNumeric(curarg) || mxIsComplex(curarg)) {
		mexErrMsgTxt(
				"register_logdetpenalty_composition: unewc must contain real arrays");
	}
	classID = mxGetClassID(curarg);

	/* Check the data type consistency between prhs[0] & prhs[1] */
	curarg = prhs[0];
	if (!mxIsEmpty(curarg)) {
		check_cell_args(curarg, d, classID, "uoldc");
	}
	check_cell_args(prhs[1], d, classID, "unewc");

	/* Call the inner wrapper mexWrapper for single or double class inputs */
	if (classID == mxSINGLE_CLASS) {
		mexWrapper<float>(nlhs, plhs, nrhs, prhs, mxSINGLE_CLASS);
	}
	else if (classID == mxDOUBLE_CLASS) {
		mexWrapper<double>(nlhs, plhs, nrhs, prhs, mxDOUBLE_CLASS);
	}
	else {
		mexErrMsgTxt(
				"register_logdetpenalty_composition: supports only single and double");
	}
}

/*
 * This is the MATLAB "inner" wrapper
 */

template <class T>
void mexWrapper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],
		mxClassID classID) {
	
	int n_dims, dimIndex, n, N;
	int szu[MAX_DIMS + 1], sz[MAX_DIMS + 1];
	const T *unew[MAX_DIMS], *uold[MAX_DIMS];
	T *ucomp[MAX_DIMS], *ucompgrad[MAX_DIMS];
	T *LDetJ, *valp, *grad[MAX_DIMS];
	T c;
	const mxArray *curarg, *cellarg;
	mxArray *cellargout;
	const mwSize *sztmp;
  
	/*
	 * Parse the input arguments
	 */
	
	/* unewc (this one is guaranteed to exist, unlike uoldc) */
	curarg = prhs[1];
	n_dims = mxGetNumberOfElements(curarg);
	if (n_dims > MAX_DIMS) {
		mexErrMsgTxt("register_logdetpenalty_composition: unewc has more than the "
				"maximum number of dimensions (fix with recompile?)");
	}
	cellarg = mxGetCell(curarg, 0);
	n = mxGetNumberOfDimensions(cellarg);
	sztmp = mxGetDimensions(cellarg);
	N = mxGetNumberOfElements(cellarg);
	if (n > n_dims) {
		/* Check that any "extra" dimensions are unity */
		for (dimIndex = n_dims; dimIndex < n; ++dimIndex) {
			if (sztmp[dimIndex] != 1) {
				mexErrMsgTxt("register_logdetpenalty_composition: the shape of unewc "
						"is inconsistent with the dimensionality");
			}
		}
	}
	/* Copy the size information and pointers to each dimension of unewc */
	for (dimIndex = 0; dimIndex < n_dims; ++dimIndex) {
		szu[dimIndex] = sztmp[dimIndex];
		cellarg = mxGetCell(curarg, dimIndex);
		unew[dimIndex] = (T *) mxGetData(cellarg);
	}

	/* uoldc */
	curarg = prhs[0];
	uold[0] = NULL;
	if (!mxIsEmpty(curarg)) {
		for (dimIndex = 0; dimIndex < n_dims; ++dimIndex) {
			cellarg = mxGetCell(curarg, dimIndex);
			uold[dimIndex] = (T *) mxGetData(cellarg);
		}
	}

	/* c */
	c = 1.0;
	if (nrhs > 2) {
		c = mxGetScalar(prhs[2]);
	}

	/*
	 * Allocate space for the output arguments
	 */

	/* ucompc */
	plhs[0] = mxCreateCellMatrix(1, n_dims);
	for (dimIndex = 0; dimIndex < n_dims; ++dimIndex) {
		cellargout = mxCreateNumericArray(n_dims, szu, classID, mxREAL);
		ucomp[dimIndex] = (T*) mxGetData(cellargout);
		mxSetCell(plhs[0], dimIndex, cellargout);
	}

	/* LDetJ */
	LDetJ = NULL;
	if (nlhs > 1) {
		for (dimIndex = 0; dimIndex < n_dims; ++dimIndex) {
			if (szu[dimIndex] > 1) {
				sz[dimIndex] = szu[dimIndex] - 1;
			}
			else {
				sz[dimIndex] = 1;
			}
		}
		plhs[1] = mxCreateNumericArray(n_dims, sz, classID, mxREAL);
		LDetJ = (T*) mxGetData(plhs[1]);
	}

	/* val */
	valp = NULL;
	if (nlhs > 2) {
		plhs[2] = mxCreateNumericMatrix(1, 1, classID, mxREAL);
		valp = (T*) mxGetData(plhs[2]);
	}

	/* grad */
	grad[0] = NULL;
	if (nlhs > 3) {
		plhs[3] = mxCreateCellMatrix(1, n_dims);
		for (dimIndex = 0; dimIndex < n_dims; ++dimIndex) {
			cellargout = mxCreateNumericArray(n_dims, szu, classID, mxREAL);
			mxSetCell(plhs[3], dimIndex, cellargout);
			grad[dimIndex] = (T*) mxGetData(cellargout);
		}
		// Also allocate temporary space
		// ucompgrad
		for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
			ucompgrad[dimIndex] = (T *) mxMalloc(sizeof(T)*n_dims*n_dims*N);
			for (int i = 0; i < n_dims*n_dims*N; ++i) {
				ucompgrad[dimIndex][i] = 0.0; /* Initialize to zero */
			}
		}
	}

	/* Do the actual computation */
	if (n_dims == 1) {
		regldpc_work<T, 1>(szu, uold, unew, c, ucomp, LDetJ, valp, grad, ucompgrad);
	}
	else if (n_dims == 2) {
		regldpc_work<T, 2>(szu, uold, unew, c, ucomp, LDetJ, valp, grad, ucompgrad);
	}
	else if (n_dims == 3) {
		regldpc_work<T, 3>(szu, uold, unew, c, ucomp, LDetJ, valp, grad, ucompgrad);
	}

	// Free temporary space, if applicable
	if (grad[0] != NULL) {
		for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
			mxFree(ucompgrad[dimIndex]);
	}

}

/* Do the actual computation */
template<class T, int n_dims>
	void regldpc_work(const int szu[], const T *uold[], const T *unew[], T c,
		T *ucomp[], T *LDetJ, T *valp, T *grad[], T *ucompgrad[]) {

	int dimIdx, dimIdx2;
	T x1[MAX_DIMS], x2[MAX_DIMS], xgrad[MAX_DIMS * MAX_DIMS];
	T *ucgP[MAX_DIMS];
	T *xp;

	bool calc_val = (valp != NULL);
	bool calc_grad = (grad[0] != NULL);

	/* Set up the quadratic interpolation of u-grids */
	Qinterp::qinterp<T, n_dims, Qinterp::reflect> q(szu);
	/* Initialize the pixel iterator*/
	pixIterator pIu(szu, n_dims, false);
	for (dimIdx = 0; dimIdx < n_dims; dimIdx++) {
	    ucgP[dimIdx] = ucompgrad[dimIdx];
	}

	/* Do the composition */
	for (pIu.restart(); !pIu.at_end(); pIu++) {

		if (uold[0] == NULL) {
			/* uold is empty */
			for (dimIdx = 0; dimIdx < n_dims; dimIdx++) {
				ucomp[dimIdx][pIu] = unew[dimIdx][pIu];
				if (calc_grad) { /* special case: using unit matrix for gradient */
					for (dimIdx2 = 0; dimIdx2 < n_dims; dimIdx2++) {
						*(ucgP[dimIdx]) = (dimIdx2 == dimIdx);
						ucgP[dimIdx]++;
					}
				}
			}
		}
		else {
			/* uold is not empty, so we interpolate */
			for (dimIdx = 0; dimIdx < n_dims; dimIdx++) {
				x1[dimIdx] = unew[dimIdx][pIu] + pIu.coord(dimIdx);
			}
			/* only do the composition here, no gradient calculation needed  */
			if (!calc_grad) {
				q.value(x1,n_dims,uold,x2);
				for (dimIdx = 0; dimIdx < n_dims; dimIdx++) {
					ucomp[dimIdx][pIu] = x2[dimIdx] + unew[dimIdx][pIu];
				}
			}
			else {
				q.valgrad(x1,n_dims,uold,x2,xgrad);
				for (dimIdx = 0,xp = xgrad; dimIdx < n_dims; dimIdx++) {
					ucomp[dimIdx][pIu] = x2[dimIdx] + unew[dimIdx][pIu];
					// reshape the gradient data so that the first index
					// (ucompgrad[first][second]) corresponds to the derivative
					// component. Also add the derivative of the identity.
					for (dimIdx2 = 0; dimIdx2 < n_dims; dimIdx2++,xp++) {
						*(ucgP[dimIdx]) = *xp + (dimIdx2 == dimIdx);
						ucgP[dimIdx]++;
					}
				}
			}
		}

	}

	if (!calc_val) {
		return;
	}

	/* Calculate the penalty value and gradient of log(det(J)) */

	Jacobianu::jacobianu<T, n_dims> cellJ;
	cellJ.init(szu);
	*valp = 0.0;

	for (pIu.restart(); !pIu.at_end(); pIu++) {
		if (pIu.on_right_edge()) {
			continue;
		}
		/* Snip out a cell of ucomp */
		cellJ.stuff(pIu, n_dims, ucomp);
		/* Calculate log(det(J)) & gradient of this cell */
		cellJ.calcLDetJ(pIu, n_dims, c, grad, ucompgrad);
		/* Get the contribution to the penalty value */
		*LDetJ = cellJ.getCellLDetJ();
		LDetJ++;
		*valp += cellJ.getCellLDetJ();
	}
	cellJ.~jacobianu();

}



#ifdef MAIN
int main(int argc,char *argv[])
{
	char *filein,*fileout;
	const int n_inputs = 3;
	const int n_outputs = 4;
	mxArray *input[n_inputs];
	mxArray *output[n_outputs];
	const char *input_names[] = {
			"uold",
			"unew",
			"c"
	};
	const char *output_names[] = {
			"ucomp",
			"detJ",
			"val",
			"grad"
	};

	if (argc < 3) {
		printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables uold, unew, and c\n",argv[0]);
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
	printf("All done, success.\n");

	return EXIT_SUCCESS;
}
#endif

