

# compilation method 1
# system(sprintf('nvcc -I"%s/extern/include" --cuda "register_mismatch_nancorr_cuda.cu" --output-file "register_mismatch_nancorr_cuda.cpp" -arch=sm_20', matlabroot))
# mex register_mismatch_nancorr_cuda.cpp -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcufft

# using compilation method 2

CUDA_INCLUDE_PATH = /usr/local/cuda/include
CUDA_LINK_PATH = /usr/local/cuda/lib64
NVCC_FLAGS = -Xcompiler -fPIC -arch=sm_20 -Xptxas -v

MATLAB_INCLUDE_PATH = /usr/local/MATLAB/R2011b/extern/include
MATLAB_FLAGS = -lcudart -lcufft -cxx

COMMON_FILE_PATH = ./common
COMMON_FILES = $(COMMON_FILE_PATH)/cuda_common.h $(COMMON_FILE_PATH)/cuda_function.h \
			   $(COMMON_FILE_PATH)/cuda_gpu_selection.h $(COMMON_FILE_PATH)/cuda_handle_error.h \
			   $(COMMON_FILE_PATH)/cuda_kernel_planning.h $(COMMON_FILE_PATH)/cuda_time_wait.h

cuda_cleanup : cuda_cleanup.cu $(COMMON_FLIES)
	nvcc $(NVCC_FLAGS) -I$(MATLAB_INCLUDE_PATH) -c cuda_cleanup.cu 
	mex -I$(CUDA_INCLUDE_PATH) -L$(CUDA_LINK_PATH) $(MATLAB_FLAGS) cuda_cleanup.o
	rm cuda_cleanup.o

register_mismatch_nancorr_cuda : register_mismatch_nancorr_cuda.cu $(COMMON_FILES) 
	nvcc $(NVCC_FLAGS) -I$(MATLAB_INCLUDE_PATH) -c register_mismatch_nancorr_cuda.cu 
	mex -I$(CUDA_INCLUDE_PATH) -L$(CUDA_LINK_PATH) $(MATLAB_FLAGS) register_mismatch_nancorr_cuda.o
	rm register_mismatch_nancorr_cuda.o

all : 
	make cuda_cleanup 
	make register_mismatch_nancorr_cuda

clean : 
	rm *.mexa64
	



#!nvcc -c register_mismatch_nancorr_cuda.cu -Xcompiler -fPIC -I /usr/local/MATLAB/R2011a2/extern/include -arch=sm_20
#mex register_mismatch_nancorr_cuda.o -I/usr/local/cuda/include ...
#    -L/usr/local/cuda/lib64 -lcudart -lcufft -cxx


