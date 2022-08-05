% REGISTER_DETJ: compute determinant arrays
% Syntax:
%   detJ = register_detJ(J)
% where
%   J is a n-by-n-by-[sz] array
% calculates
%   detJ, an array of size [sz], where each element of detJ, say detJ(i),
%     is the determinant of the matrix J(:,:,i).
%
% Currently, the matrices appearing in the determinant must be 1-by-1,
% 2-by-2, or 3-by-3.
%
% See also: REGISTER_JACOBIAN, REGISTER_G2DETJ.

% This is implemented as a MEX file, for speed.
