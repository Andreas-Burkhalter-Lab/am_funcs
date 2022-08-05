% GRADWLLCV: Gradient of log-likelihood for gaussian kernel-density estimator
% Usage:
%  grad = gradwllcv(x,y,n,G,Pcv,mode)
% where
%   x is a d-by-N matrix of points in d-dimensional space;
%   y is a d'-by-N matrix of points in d'-dimensional space;
%   n is the multiplicity of each point;
%   G is a N-by-N matrix of gaussian kernel functions;
%   Pcv is the cross-validated density at each point (if you don't want
%     cross-validation, instead supply P);
%   mode is a string, 'full' or 'diag', depending on whether W is a
%     full matrix or a diagonal matrix
% and
%   grad is a d'-by-d matrix (full case) or d'-by-1 vector (diag case),
%     yielding the gradient of SCV with respect to the elements of W.
%
% Note that this returns the gradient of
%   -n*log(Pcv)
% and does not include the -N*log(abs(det(W*W'))/2 term. The derivative
% of that term is simply
%    -N*inv(W*W')*W
%
% This is written in C for speed. The large number of arguments are
% passed simply to avoid re-calculation.
