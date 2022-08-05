function L1 = calcL1d1(dx)
% CALCL1D1: calculate the matrix in discrete action for d=1.
% Syntax:
%   L1 = calcL1d1(dx)
% where
%   dx is the vector of separations, dx_i = x_{i+1}-x{i}.
% and
%   L1 is a tridiagonal matrix containing the kinetic energy & volume
%     operators.  It's the inverse of the matrix
%              G_{ij} = exp(-|x_i - x_j|)/2.
%
% See also: LLCV1, CALCVD1.

% Copyright 2006 by Timothy E. Holy

dx = dx(:);
sdx = sinh(dx);
od = -1./sdx;	% The off-diagonal elements
dw = cosh(dx)./sdx;
dw(isnan(dw)) = 1;  % When points are widely separated, cosh(dx)/sinh(dx)=1
d = [1;dw] + [dw;1];
% The sparse matrix formulation has too much overhead to justify using
% it---it can dominate the calculation---and so use a dedicated
% tridiagonal solver instead. Just return the diagonals
L1 = {od,d,od};
%M = length(dx)+1;
%odi = 1:(M-1);
%di = 1:M;
%L1 = sparse([odi,odi+1,di],[odi+1,odi,di],[od;od;d],M,M);
%return
