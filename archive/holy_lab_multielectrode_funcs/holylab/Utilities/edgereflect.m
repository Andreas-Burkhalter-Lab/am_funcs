function xo = edgereflect(sz,x)
% EDGEREFLECT: coordinates of arrays with mirror-reflection symmetry
% In interpolation and extrapolation one sometimes wants to access array
% values that lie beyond the edge. This function implements one mechanism
% for extrapolating, by mirror-reflecting coordinates across array
% boundaries.
%
% Syntax:
%   xo = edgereflect(sz,x)
% where
%   sz is a 1-by-n_dims vector giving the size of the array
%   x is an n_pts-by-n_dims matrix specifying the coordinates of points,
%     which may extend beyond the edge of the array
% and
%   xo is an n_pts-by-n_dims matrix with coordinates that are within the
%     bounds of the array.

% Copyright 2011 by Timothy E. Holy

  xo = bsxfun(@mod,x-1,2*sz)+1;
  for dimIndex = 1:length(sz)
    lookup = [1:sz(dimIndex) sz(dimIndex):-1:1];
    xo(:,dimIndex) = lookup(xo(:,dimIndex));
  end
  