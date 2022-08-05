%
% Mex function used to reset cuda device and release all the GPU memory resources
% 
% Syntax in matlab code:
% 	cuda_cleanup();
%
% Execute this mex function only once after completing all the cuda accelerated computations.
% Do not call this mex function repeatly in a for-loop, which 
% will significantly slow down the overall performance.
%

% By Jian Wang in Holy's Lab
