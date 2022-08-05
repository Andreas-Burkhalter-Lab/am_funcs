	
	//
	// cuda event tracking start
	//
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	
	//
	// cuda event tracking end
	//

	cudaEventRecord(stop, 0);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("time: %f \n", elapsedTime);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
