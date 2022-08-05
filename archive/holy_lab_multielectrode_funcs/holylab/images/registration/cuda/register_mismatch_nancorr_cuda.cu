// C++ header
#include <stdio.h>
#include <time.h>

// cuda header
#include <cuda.h>
#include <cufft.h>
#include <math_functions.h>
#include <cuComplex.h>

// matlab header
#include "mex.h"
#include "matrix.h"

// local header
#include "common/cuda_common.h"

//
// cuda kernel 1 declaration
//

__global__ void kernel_1( // compute thetaf, fixed, w, wf, wf2, wthetaf
		cufftDoubleReal *fixed_db,
		cufftDoubleReal *w_db,
		cufftReal *fixed, 
		cufftReal *w,
        cufftReal *wthetaf, 
        cufftReal *wf, 
        cufftReal *wf2,
        int array_length) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = blockDim.x * gridDim.x;
	
	while (tid < array_length) {
		float local_fixed = __double2float_rn(fixed_db[tid]);
		float local_w = __double2float_rn(w_db[tid]);
		int thetaf = isnan(local_fixed);
		if (thetaf == 1) local_fixed = 0.0;
		thetaf = !thetaf;
		fixed[tid] = local_fixed;
		w[tid] = local_w;
		wthetaf[tid] = local_w * __int2float_rn(thetaf);
		float local_wf = local_w * local_fixed; 
		wf[tid] = local_wf;
		wf2[tid] = local_wf * local_fixed;
		tid += offset;
	}
}

//
// cuda kernel 11 declaration
//

__global__ void kernel_11( // calculate sum(w(:)) using 1 block of 1024 threads
		cufftDoubleReal *w_db,
		int array_length,
		cufftReal *sum_w) {
		
	const uint block_sum_length = 1024;
	__shared__ cufftReal block_sum[block_sum_length];
	block_sum[threadIdx.x] = 0.0;
	
	int tid = threadIdx.x;
	int offset = blockDim.x;
	while (tid < array_length) {
		block_sum[threadIdx.x] += __double2float_rn(w_db[tid]);
		tid += offset;
	}
	uint i = block_sum_length / 2;
	while (i > blockDim.x) {
		i /= 2;
	}
	while (i != 0) {
		__syncthreads();
		if (threadIdx.x < i)
			block_sum[threadIdx.x] += block_sum[threadIdx.x + i];
		i /= 2;
	}
	if (threadIdx.x == 0)
		sum_w[0] = block_sum[0];
}


//
// cuda kernel 2 declaration
//

__global__ void kernel_2(
		cufftDoubleReal *moving_db,
		cufftReal *moving, 
		cufftReal *thetam,
        cufftReal *m2, 
        int array_length) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x; // compute moving, thetam, nanflag, m2
	int offset = blockDim.x * gridDim.x;
	while (tid < array_length) {
		float local_moving = __double2float_rn(moving_db[tid]);
		int temp = isnan(local_moving);
		if (temp == 1)
			local_moving = 0.0;
		temp = !temp;
		thetam[tid] = __int2float_rn(temp);
		moving[tid] = local_moving;
		m2[tid] = local_moving * local_moving;
		tid += offset;
	}
}

//
// cuda kernel 3 declaration
//

__global__ void kernel_3(
		cufftComplex *wf_fft, 
		cufftComplex *m_fft,
        cufftComplex *wthetaf_fft, 
        cufftComplex *m2_fft,
        cufftComplex *wf2_fft,
        cufftComplex *thetam_fft,
        cufftComplex *numerator_fft,
        cufftComplex *denominator_fft,
        int array_length,
        cufftReal *sum_w) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x; // compute numerator and denominator before ifftn
	int offset = blockDim.x * gridDim.x;
	while (tid < array_length) {
		cufftComplex c1 = make_cuFloatComplex(-2.0, 0.0);
		c1 = cuCmulf(c1, cuConjf(wf_fft[tid]));
		c1 = cuCmulf(c1, m_fft[tid]);	
		cufftComplex c2 = cuCmulf(cuConjf(wthetaf_fft[tid]), m2_fft[tid]);
		c1 = cuCaddf(c1, c2);
		c2 = cuCmulf(cuConjf(wf2_fft[tid]), thetam_fft[tid]);
		numerator_fft[tid] = cuCaddf(c1, c2);	
		
		c1 = cuCmulf(cuConjf(wthetaf_fft[tid]), thetam_fft[tid]);
		c2 = make_cuComplex((1.0 / sum_w[0]), 0.0);
		denominator_fft[tid] = cuCmulf(c1, c2);
		
		tid += offset;
	}
}

//
// cuda kernel_fftshift declaration
//

__global__ void kernel_fftshift(
		cufftDoubleReal *numerator_db,
		cufftDoubleReal *denominator_db,
		cufftReal *numerator,
		cufftReal *denominator,
		int array_length) {
	
	dim3 offset; // array size along each dimension
	offset.x = gridDim.x * blockDim.x;
	offset.y = gridDim.y * blockDim.y;
	offset.z = gridDim.z * blockDim.z;
	
	dim3 from; // array index: from
	from.x = threadIdx.x + blockDim.x * blockIdx.x;
	from.y = threadIdx.y + blockDim.y * blockIdx.y;
	from.z = threadIdx.z + blockDim.z * blockIdx.z;
	
	dim3 to; // array index: to
	(from.x > (offset.x / 2 - 1)) ? (to.x = from.x - offset.x / 2) : (to.x = from.x + offset.x / 2);
	(from.y > (offset.y / 2 - 1)) ? (to.y = from.y - offset.y / 2) : (to.y = from.y + offset.y / 2);
	(from.z > (offset.z / 2 - 1)) ? (to.z = from.z - offset.z / 2) : (to.z = from.z + offset.z / 2);
	
	int from_tid = static_cast<int> (from.z + from.y * offset.z + from.x * offset.y * offset.z);
	int to_tid = static_cast<int> (to.z + to.y * offset.z + to.x * offset.y * offset.z);
	
	// devide numerator and denominator after ifftn with array_length
	cufftReal array_length_flt = static_cast<float>(array_length);
	
	numerator_db[to_tid] = static_cast<double>(numerator[from_tid] / array_length_flt);
	denominator_db[to_tid] = static_cast<double>(denominator[from_tid] / array_length_flt);
}



//
// mexFunction
//

void mexFunction(int nlhs, mxArray *plhs[], 
							int nrhs, const mxArray *prhs[]) {

	//
	// check and initialize the input & output arguments
	//

	if (nrhs != 3) // check for proper number of input arguments 
		mexErrMsgTxt("Three input arguments are required.");
	if (nlhs != 2) // check for proper number of output arguments
		mexErrMsgTxt("Two output arguments are required.");

	for (int i = 0; i < nrhs; ++i) {
		if (!mxIsDouble(prhs[i])) // make sure all the input arguments are double
			mexErrMsgTxt("All input arguments must be double.");
	}

	double *h_fixed_db = mxGetPr(prhs[0]); // h_ prefix implies host variable
															        // _db suffix implies double precission
																	// other variables are single precission by default
	mwSize h_fixed_dim = mxGetNumberOfDimensions(prhs[0]);
	const mwSize *h_fixed_size = mxGetDimensions(prhs[0]);

	double *h_moving_db = mxGetPr(prhs[1]);
	mwSize h_moving_dim = mxGetNumberOfDimensions(prhs[1]);
	const mwSize *h_moving_size = mxGetDimensions(prhs[1]);

	double *h_w_db = mxGetPr(prhs[2]);
	mwSize h_w_dim = mxGetNumberOfDimensions(prhs[2]);
	const mwSize *h_w_size = mxGetDimensions(prhs[2]);
	
	for (int i = 0; i < h_fixed_dim; ++i) {
		if ((h_fixed_size[i] != h_moving_size[i])
		        || (h_moving_size[i] != h_w_size[i])
		        || (h_w_size[i] != h_fixed_size[i]))
			mexErrMsgTxt("The input arguments are not of the same size.");
	}

	mwSize plhs_dim = h_fixed_dim; // create output arguments
	const mwSize *plhs_size = h_fixed_size;
	plhs[0] = mxCreateNumericArray(plhs_dim, plhs_size, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(plhs_dim, plhs_size, mxDOUBLE_CLASS, mxREAL);
	double *h_numerator_db = mxGetPr(plhs[0]);
	double *h_denominator_db = mxGetPr(plhs[1]);
	
	//
	// define variables related to general array size
	//
	
	int array_dim = static_cast<int>(h_fixed_dim);
	int *array_size = (int *) malloc(array_dim * sizeof(int)); // univeral array size
	int array_length = 1;
	for (int i = 0; i < array_dim; ++i) {
		// switch the array size as matlab & c are using different array storage orders
		array_size[i] = static_cast<int>(h_fixed_size[array_dim - 1 - i]);
		array_length *= array_size[i];
	}
	size_t rBytes = static_cast<size_t>(array_length) * sizeof(cufftReal); // byte size of real array
	size_t dBytes = rBytes * static_cast<size_t>(2); // byte size of double precission array
	
	//
	// define variables related to "fixed" & "w"
	//
	
	cufftDoubleReal *d_fixed_db, *d_w_db; // double precission variables
	cudaMalloc((void **) &d_fixed_db, dBytes);
	cudaMalloc((void **) &d_w_db, dBytes);	
	
	cufftReal *d_fixed, *d_w, *d_wthetaf, *d_wf, *d_wf2;	// single precission variables
	cudaMalloc((void **) &d_fixed, rBytes);
	cudaMalloc((void **) &d_w, rBytes);
	cudaMalloc((void **) &d_wthetaf, rBytes);
	cudaMalloc((void **) &d_wf, rBytes);
	cudaMalloc((void **) &d_wf2, rBytes);

	//
	// define variables related to "moving"
	//
	
	cufftDoubleReal *d_moving_db; // double precission variables
	cudaMalloc((void **) &d_moving_db, dBytes);

	cufftReal *d_moving, *d_thetam, *d_m2; // single precission variables
	cudaMalloc((void **) &d_moving, rBytes);
	cudaMalloc((void **) &d_thetam, rBytes);
	cudaMalloc((void **) &d_m2, rBytes);
	
	//
	// copy double-precission variables from Host to Device
	//
	
	cudaMemcpy(d_fixed_db, h_fixed_db, dBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_w_db, h_w_db, dBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_moving_db, h_moving_db, dBytes, cudaMemcpyHostToDevice);
	
	//
	// pre-computation related to "fixed" & "w", "moving" before fftn
	//
	
	dim3 blocks, threads, size(array_length, 1, 1);
	cuda_kernel_planning(blocks, threads, size); // optimize these, block size calculator
	
	cufftReal *d_sum_w;
	cudaMalloc((void **) &d_sum_w, sizeof(cufftReal));
	
	kernel_11<<<1, 1024>>> (
			d_w_db,
			array_length,
			d_sum_w);

	kernel_1<<<blocks, threads>>> (
			d_fixed_db,
			d_w_db,
			d_fixed,
			d_w,
			d_wthetaf,
			d_wf,
			d_wf2,
			array_length);

	kernel_2<<<blocks, threads>>> (
			d_moving_db,
			d_moving,
			d_thetam,
			d_m2,
			array_length);
	
	
	cudaThreadSynchronize(); // make sure all kernels are completed
	cuda_handle_error();
	
	cudaFree(d_fixed_db);
	cudaFree(d_w_db);
	cudaFree(d_moving_db);
	cudaFree(d_fixed);
	cudaFree(d_w);

	//
	// compute the individual terms for the quadratic
	//

	int array_length_cplx = 1; // complex array length which is less than array_length
	for (int i = 0; i < (array_dim-1); ++i)
		array_length_cplx *= array_size[i];
	array_length_cplx *= (array_size[array_dim-1] / 2 + 1);	
	size_t cBytes = static_cast<size_t>(array_length_cplx) * sizeof(cufftComplex);

	cufftComplex *d_wthetaf_fft, *d_m2_fft, *d_wf_fft, *d_m_fft;
	cudaMalloc((void**) &d_wthetaf_fft, cBytes);
	cudaMalloc((void**) &d_m2_fft, cBytes);
	cudaMalloc((void**) &d_wf_fft, cBytes);
	cudaMalloc((void**) &d_m_fft, cBytes);

	cufftComplex *d_wf2_fft, *d_thetam_fft;
	cudaMalloc((void **) &d_wf2_fft, cBytes);
	cudaMalloc((void **) &d_thetam_fft, cBytes);

	int istride = 1;
	int idist = 0;
	int ostride = 1;
	int odist = 0;
	int batch = 1;

	cufftHandle plan; // create cufft plan
	cufftPlanMany(
			&plan, array_dim, array_size, 
			NULL, istride, idist, 
			NULL, ostride, odist, 
			CUFFT_R2C, batch);
	
	cufftExecR2C(plan, d_wthetaf, d_wthetaf_fft);
	cufftExecR2C(plan, d_m2, d_m2_fft);
	cufftExecR2C(plan, d_wf, d_wf_fft);
	cufftExecR2C(plan, d_moving, d_m_fft);
	
	cufftExecR2C(plan, d_wf2, d_wf2_fft);
	cufftExecR2C(plan, d_thetam, d_thetam_fft);

	cudaThreadSynchronize(); // make sure all FFTN are completed
	cuda_handle_error();
	cufftDestroy(plan);
	
	cudaFree(d_wthetaf);
	cudaFree(d_m2);
	cudaFree(d_wf);
	cudaFree(d_moving);
	
	cudaFree(d_wf2);
	cudaFree(d_thetam);
	
	//
	// assemble the quadratic and compute ifftn terms
	//

	cufftComplex *d_numerator_fft, *d_denominator_fft;
	cudaMalloc((void **) &d_numerator_fft, cBytes);
	cudaMalloc((void **) &d_denominator_fft, cBytes);
	
	kernel_3<<<blocks, threads>>> (
			d_wf_fft, 
			d_m_fft,
	        d_wthetaf_fft, 
	        d_m2_fft,
	        d_wf2_fft,
	        d_thetam_fft,
	        d_numerator_fft,
	        d_denominator_fft,
	        array_length_cplx,
	        d_sum_w);
	
	cudaThreadSynchronize(); // make sure kernel_3 is completed
	cuda_handle_error();
	
	cudaFree(d_wthetaf_fft);
	cudaFree(d_m2_fft);
	cudaFree(d_wf_fft);
	cudaFree(d_m_fft);
	cudaFree(d_wf2_fft);
	cudaFree(d_thetam_fft);
	cudaFree(d_sum_w);
	
	cufftReal *d_numerator, *d_denominator;
	cudaMalloc((void **) &d_numerator, rBytes);
	cudaMalloc((void **) &d_denominator, rBytes);
	
	cufftPlanMany(&plan, 
			array_dim, array_size, 
			NULL, istride, idist, 
			NULL, ostride, odist, 
			CUFFT_C2R, batch);;
	
	cufftExecC2R(plan, d_numerator_fft, d_numerator);
	cufftExecC2R(plan, d_denominator_fft, d_denominator);
	
	cudaThreadSynchronize(); // make sure all IFFTN are completed
	cuda_handle_error();
	cufftDestroy(plan);
	
	cudaFree(d_numerator_fft);
	cudaFree(d_denominator_fft);
	
	//
	// norminaze denominator & numerator
	//
	
	cufftDoubleReal *d_numerator_db, *d_denominator_db;
	cudaMalloc((void **) &d_numerator_db, dBytes);
	cudaMalloc((void **) &d_denominator_db, dBytes);
	
	size.x = static_cast<uint> (array_size[0]);
	size.y = static_cast<uint> (array_size[1]);
	size.z = static_cast<uint> (array_size[2]);
	cuda_kernel_planning2(blocks, threads, size);
	
	kernel_fftshift<<<blocks, threads>>> (
				d_numerator_db,
				d_denominator_db,
				d_numerator,
				d_denominator,
				array_length);
	
	cudaThreadSynchronize();
	cudaMemcpy(h_numerator_db, d_numerator_db, dBytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_denominator_db, d_denominator_db, dBytes, cudaMemcpyDeviceToHost);
	cuda_handle_error();

	cudaFree(d_numerator);
	cudaFree(d_denominator);
	cudaFree(d_numerator_db);
	cudaFree(d_denominator_db);
	
	//
	// release allocated cpu memory
	//
	
	free(array_size);
}



