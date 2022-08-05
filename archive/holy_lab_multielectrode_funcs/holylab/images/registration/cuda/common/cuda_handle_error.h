

#ifndef __CUDA_HANDLE_ERROR_H__
#define __CUDA_HANDLE_ERROR_H__

// system header
#include <stdio.h>

// cuda header
#include <cuda.h>
#include <cufft.h>

// cuda runtime API error handling
#define cuda_handle_error() { \
	cudaError_t cuda_error = cudaGetLastError(); \
	if (cuda_error != cudaSuccess) \
		printf("Error '%s' in %s: line %d \n", cudaGetErrorString(cuda_error), __FILE__, __LINE__); \
}

#endif  // __CUDA_HANDLE_ERROR_H__
