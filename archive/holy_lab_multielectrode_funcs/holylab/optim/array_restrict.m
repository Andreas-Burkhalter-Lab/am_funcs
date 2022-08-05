% ARRAY_RESTRICT: generate a coarse-grid version of an array (useful for multigrid)
%
% This function is the "inverse" of ARRAY_PROLONG, which uses
% interpolation to generate a fine-scale version of an array. This
% restriction function satisfies the Galerkin condition, meaning that
% the restriction operator is (up to a constant factor) the transpose of
% the prolongation operator (the constant is chosen to preserve the mean
% value).  As a consequence, this function properly de-aliases (smooths)
% the data in the course of restriction.
%
% Syntax:
%   Ar = array_restrict(A)
% where A is the fine-grid array, and Ar is the coarse-grid version of
% the same array.
%
%   Ar = array_restrict(A,dimFlag)
% where dimFlag is a logical vector of length ndims(A), where dimensions
% set to false are not coarsened.
%
%   Ar = array_restrict(A,dimFlag,n_threads)
% where n_threads is an integer specifying the number of threads to use
% (NOTE: currently disabled, the single-threaded version is better).
%
% A special feature of this function is that it can restrict arrays of
% any size, including along dimensions that have an even size.  See
% ARRAY_PROLONG for details.
%
% See also: ARRAY_PROLONG.

% Copyright 2009 by Timothy E. Holy
% Implemented as a MEX file for speed
