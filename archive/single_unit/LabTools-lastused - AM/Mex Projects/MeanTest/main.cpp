#include <mex.h>
#include <matrix.h>


void mexFunction(
	int nlhs,              // Number of left hand side (output) arguments
	mxArray *plhs[],       // Array of left hand side arguments
	int nrhs,              // Number of right hand side (input) arguments
	const mxArray *prhs[]  // Array of right hand side arguments
)
{
	int ok;
	double sizeOfVector;
	mxArray *x[1];


	if (nlhs != 1) {
		mexErrMsgTxt("Invalid # of return values.\n");
		return;
	}

	if (nrhs != 1) {
		mexErrMsgTxt("Invalid # of args.\n");
		return;
	}


	// Copy the vector over into another mxArray.
	x[0] = mxDuplicateArray(prhs[0]);

	mexCallMATLAB(nlhs, plhs, 1, x, "mean");
}