function V = calcVd1(dx)
% CALCVD1: calculate the volume matrix for the continuous d=1 theory
% Syntax:
%   V = calcVd1(dx)
% where
%   dx is the vector of separations, dx_i = x_{i+1}-x{i}.
%
% See also: CALCL1D1.

% Copyright 2006 by Timothy E. Holy

dx = dx(:);
M = length(dx)+1;
sdx = sinh(dx);
cdx = cosh(dx);
sdx2 = sdx.^2;
od = -1./sdx + dx.*cdx./sdx2;
d = cdx./sdx - dx./sdx2;
d(isnan(d)) = 1;  % For widely-separated points, cdx & sdx are Inf
od(isnan(od)) = 0;
% The sparse formalism takes too much overhead, use dedicated solver.
V = {od/2,([d;1]+[1;d])/2,od/2};
%B = [[od;0],[d;1]+[1;d],[0;od]]/2;
%V = spdiags(B,[-1 0 1],M,M);
