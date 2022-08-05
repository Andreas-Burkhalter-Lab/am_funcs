%
% register_mismatch_nancorr_cuda : CUDA version of register_mismatch_nancorr()
% 
% This code uses a similar syntax as register_mismatch_nancorr() and will be
% called directly by register_block_mismatch_cuda() by far.
%
% After finishing image registration using CUDA code,
% call cuda_cleanup() to release & reset CUDA GPU

% By Jian Wang in Holy Lab
