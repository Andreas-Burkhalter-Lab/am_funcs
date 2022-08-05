/*
 * register_barrierpenalty_composition.cpp : calculate the edge barrier penalty
 * inside of register_block_penalty.m
 *
 * Copyright 2012 Jian Wang
 */

#include "mex.h"
#include "matrix.h"
#include <cmath> 						/* for log */
#include <limits> 					/* for NaN */

template<class T>
void barrierPenaltyWork(T *unew, int unewNumElem, int mismatchNumDims,
			const int *mismatchDims, T *BPsum, T *BPgrad);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Check the input & output argument numbers */
	if (nrhs != 2) {
		mexErrMsgTxt("register_barrierpenalty_composition: requires two inputs");
	}
	if (nlhs < 1 || nlhs > 2) {
		mexErrMsgTxt("register_barrierpenalty_composition: must have 1 or 2 outputs");
	}

	/* Parse the input arguments */
	/* prhs[0] */
	const mxArray *curarg = prhs[0];
	if (!mxIsNumeric(curarg) || mxIsComplex(curarg)) {
	    mexErrMsgTxt("register_barrierpenalty_composition: unew must be a real array");
	}
	mxClassID classID = mxGetClassID(curarg);
	int unewNumDims = mxGetNumberOfDimensions(curarg);
	const int *unewDims = mxGetDimensions(curarg);
	int unewNumElem = mxGetNumberOfElements(curarg);
	void *unew = mxGetData(curarg);

	/* prhs[1] */
	curarg = prhs[1];
	if (!mxIsNumeric(curarg) || mxIsComplex(curarg)) {
		mexErrMsgTxt("register_barrierpenalty_composition: mismatch must be a real array");
	}
	int mismatchNumDims = mxGetNumberOfDimensions(curarg);
	const int *mismatchDims = mxGetDimensions(curarg);

	/* Allocate space for the output arguments */
	/* plhs[0] */
	plhs[0] = mxCreateNumericMatrix(1,1,classID,mxREAL);
	void *BPsum = mxGetData(plhs[0]);
	
	/* plhs[1] */
	void *BPgrad = NULL;
	if (nlhs == 2) {
		plhs[1] = mxCreateNumericArray(unewNumDims,unewDims,classID,mxREAL);
		BPgrad = mxGetData(plhs[1]);
	}

	/* Do the actual computation */
	if (classID == mxSINGLE_CLASS) {
		barrierPenaltyWork<float>((float *) unew, unewNumElem, mismatchNumDims,
				mismatchDims, (float *) BPsum, (float *) BPgrad);
	}
	else if (classID == mxDOUBLE_CLASS) {
		barrierPenaltyWork<double>((double *) unew, unewNumElem, mismatchNumDims,
				mismatchDims, (double *) BPsum, (double *) BPgrad);
	}
}


template<class T>
void barrierPenaltyWork(T *unew, int unewNumElem, int mismatchNumDims,
			const int *mismatchDims, T *BPsum, T *BPgrad) {

	bool calc_grad = (BPgrad != NULL);

	for (int i = 0; i < unewNumElem; ++i) {
		T gridSize = mismatchDims[i];
		T center = ceil(gridSize/2.0);
		T coord = center + unew[i];
		T gridGap[2] = {1.5, 0.5};
		if (coord < gridGap[0] || coord > (gridSize - gridGap[1])) {
			*BPsum = std::numeric_limits<T>::quiet_NaN();
		}
		else {
			*BPsum = *BPsum + log(coord - gridGap[0]) + log(gridSize - gridGap[1] - coord);
			if (calc_grad) {
				BPgrad[i] = 1.0/(coord - gridGap[0]) - 1.0/(gridSize - gridGap[1] - coord);
			}
		}
	}

}

