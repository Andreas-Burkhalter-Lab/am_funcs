
#ifndef __CUDA_KERNEL_PLANNING_H__
#define __CUDA_KERNEL_PLANNING_H__

// system header
#include <stdio.h>

// cuda header
#include <cuda.h>

// local header
#include "cuda_gpu_selection.h"

#define block_size_x 256

//
// cuda kernel lauch planning function (beta version, 1D array)
//

void cuda_kernel_planning(
		dim3 &blocks, 
		dim3 &threads,
		dim3 &length) {
	
	if ((length.y == 1) && (length.z == 1)) {
		threads.x = block_size_x; // 8 resident blocks per multiprocessor
		threads.y = 1;
		threads.z = 1;
		
		blocks.x = (length.x + threads.x - 1) / threads.x;
		if (blocks.x > gpu_prop_maxGridSize_0) {
			//printf("1-D block is not enough. \n");
			blocks.x = gpu_prop_maxGridSize_0;
		}
		blocks.y = 1;
		blocks.z = 1;
		
		//printf("inside of kernel planning: %d %d %d \n", length.x, blocks.x, threads.x);
	}
	
}

void cuda_kernel_planning2(
		dim3 &blocks, 
		dim3 &threads,
		dim3 &arrayDims) {
	
	// each dim of image block are ceiled up on a 32 base
	// 3-d thread block: 1 x 16 x 16 (x, y, z)
	threads.x = 1;
	threads.y = 16;
	threads.z = 16;

	// 3-d block grid:
	blocks.x = (arrayDims.x + threads.x - 1) / threads.x;
	blocks.y = (arrayDims.y + threads.y - 1) / threads.y;
	blocks.z = (arrayDims.z + threads.z - 1) / threads.z;
}

#endif  // __CUDA_KERNEL_PLANNING_H__
