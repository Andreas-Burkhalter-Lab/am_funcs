function [P,PCV,dPdX,dPdz] = pderivatives(X,n)
% PDERIVATIVES: calculate P, PCV, and useful derivatives
% Syntax:
%   [P,PCV,dPdX,dPdz] = pderivatives(X,n)
% where
%   X are the projected data points, a M-by-dp (dp = # projection
%     dimensions) matrix;
%   n contains the multiplicity of each data point, a 1-by-M vector;
% and all of the following are of size M-by-dp (so each column corresponds
% to a projection direction):
%   P is the kernel estimate, using exponential kernels;
%   PCV is the cross-validated kernel estimate;
%   dPdX is the slope of P at each data point;
%   dPdz is the derivative of P (and PCV) with respect to the global scale,
%     i.e., dPdz = (d/dz) P(z*X), evaluated at z=1.

% Copyright 2006 by Timothy E. Holy

% Possible further options?
%   sorted
%   c, s, & cholL1 precalculated?

if (nargin < 4)
  options = struct;
end

N = sum(n);
[M,dp] = size(X);
if (M ~= length(n))
  error('The sizes of X and n do not match');
end

% Project the data points & sort the projections
[X,sort_order] = sort(X);
n = n(sort_order);
dX = diff(X);
npN = n/N;
if (dp == 1)
  npN = npN';
end

c = cosh(dX);
s = sinh(dX);
index_inf = find(c == Inf);   % In cases where points are widely separated
d0 = c./s;
d0(index_inf) = 1;

% Calculate P, PCV
% Here I use a dedicated tridiagonal solver I wrote. It makes it
% faster than using Matlab's sparse matrices. It also makes it quite a bit
% more convenient to do multiple projections at once
v0 = zeros(1,dp);
v1 = ones(1,dp);
od = -1./s;
d1 = [v1;d0];
d2 = [d0;v1];
dd = d1 + d2;
P = tridiag_inv(od,dd,od,npN); % P = L\npN
if any(isnan(P))
  error('P has nans');
end
%PCV = P - npN/2;
PCV = crossval1d(P',dX)';

% Slope of P
% Since the slope is not continuous across the data points, we
% choose this as the average of the slopes on either side of the data
% points.
dPdX = ([P(2:end,:);v0]./[s;v1] - P.*d2  + P.*d1 - [v0;P(1:end-1,:)]./[v1;s])/2;

% Derivative of P (& PCV) with respect to global scale
% First compute \partial L/\partial z, where z is global scale, and
% multiply it by P
divs2 = dX./s.^2;
cdivs2 = c.*divs2;
cdivs2(index_inf) = 0;
dLdzP = [v0;cdivs2].*[v0;P(1:end-1,:)] - ([divs2;v0] + [v0;divs2]).*P + ...
        [cdivs2;v0].*[P(2:end,:);v0];
% Now finish it by dividing by L
dPdz = -tridiag_inv(od,dd,od,dLdzP);
if any(isnan(dPdz))
  error('dPdz has nans');
end

% Finally, invert the sort permutation
[sso,sort_invert] = sort(sort_order);
for i = 1:dp
  P(:,i) = P(sort_invert(:,i),i);
  PCV(:,i) = PCV(sort_invert(:,i),i);
  dPdX(:,i) = dPdX(sort_invert(:,i),i);
  dPdz(:,i) = dPdz(sort_invert(:,i),i);
end