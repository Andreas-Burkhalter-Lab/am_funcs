function v = tridiag_inv(dm1,d,dp1,r)
% TRIDIAG_INV: solve tridiagonal systems
% For the linear equation
%     Ax = r,
% where A is an n-by-n tridiagonal matrix, solve for x.
%
% Syntax:
%   x = tridiag_inv(dm1,d,dp1,r)
% where
%   d is a 1-by-n vector containing the diagonal of A;
%   dm1 is a 1-by-(n-1) vector containing the lower off-diagonal of A;
%   dp1 is a 1-by-(n-1) vector containing the upper off-diagonal of A;
%   r is a 1-by-n vector containing the right hand side;
% and
%   x is a 1-by-n vector containing the solution to the linear equation.
%
% Alternatively, all can be matrices. In that case, each corresponding
% set of columns is an independent tridiagonal system.
%
% This is a MEX file for speed.
%
% See also: TRIDIAG_MULT.