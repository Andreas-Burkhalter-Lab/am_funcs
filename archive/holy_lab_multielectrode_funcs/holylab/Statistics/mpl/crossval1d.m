function ucv = crossval1d(u,dX)
% CROSSVAL1D: cross-validated density/amplitude in d=1
% Cross-validation calculates the density (Pcv) or amplitude (psicv) at
% the current point using all points but the current one. This function
% does some of the arithmetic analytically to avoid roundoff error. The
% expression is
%    ucv = (1 - L/2)*u
% where u is either psi or P. This function calculates 1-L/2 in a way that
% is relatively insensitive to roundoff error.
%
% Syntax:
%   psicv = crossval1d(psi,dX)
%   Pcv = crossval1d(P,dX)
% where
%   psi/P is the amplitude/density (depending on whether you're doing
%     MPL or kernel estimation)
%   dX is diff(X), where X contains the projected 1-d coordinate, assumed
%     to be sorted
% and
%   psicv/Pcv is the cross-validated amplitude/density
  
% Copyright 2006 by Timothy E. Holy
  
  %dX = diff(X);
  edX = exp(-dX);
  sdX = sinh(dX);
  eos = edX./(2*sdX);
  s2 = 1./(2*sdX);
  ucv = [s2.*u(2:end), 0] + [0, s2.*u(1:end-1)] - [eos.*u(1:end-1), 0] ...
        - [0, eos.*u(2:end)];
  