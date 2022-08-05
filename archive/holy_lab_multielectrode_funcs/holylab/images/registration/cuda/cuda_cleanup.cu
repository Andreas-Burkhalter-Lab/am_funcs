// C++ header
#include <stdio.h>

// cuda header
#include <cuda.h>

// matlab header
#include "mex.h"
#include "matrix.h"

// local header
#include "common/cuda_common.h"

//
// mexFunction: entrance point
//

void mexFunction(int nlhs, mxArray *plhs[], 
							int nrhs, const mxArray *prhs[]) {
	
	//
	// release allocated cuda gpu device
	// execute once after finishing cuda part calculation
	//
	
	cudaDeviceReset();
	
}

