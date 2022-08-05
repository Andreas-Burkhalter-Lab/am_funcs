function r = tridiag_mult(dm1,d,dp1,x)
% TRIDIAG_MULT: multiply by a tridiagonal matrix
% This calculates Ax, where x is a vector and A a tridiagonal matrix.
%
% Syntax:
%   r = tridiag_mult(dm1,d,dp1,x)
% where
%   d is a 1-by-n vector containing the diagonal of A;
%   dm1 is a 1-by-(n-1) vector containing the lower off-diagonal of A;
%   dp1 is a 1-by-(n-1) vector containing the upper off-diagonal of A;
%   x is a 1-by-n vector containing the vector to be multiplied
% and
%   r is a 1-by-n vector containing the right hand side;
%
% Alternatively, all can be matrices. In that case, each corresponding
% set of columns is an independent tridiagonal system.
%
% See also: TRIDIAG_INV.
  
% Copyright 2006 by Timothy E. Holy
  
% Note that this could be made faster with a MEX file, since the slow
% step will be the allocation of memory for the intermediates.
  
  n_problems = size(x,2);
  z = zeros(1,n_problems);
  r = [z; dm1 .* x(1:end-1,:)] + d .* x + [dp1 .* x(2:end,:); z];
  