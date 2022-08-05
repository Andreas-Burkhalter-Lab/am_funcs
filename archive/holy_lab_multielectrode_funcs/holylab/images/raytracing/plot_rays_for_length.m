function hline = plot_rays_for_length(r,l)
% plot_rays_for_length: plot rays for a certain propagation distance
% Syntax:
%   hline = plot_rays_for_length(r,l)
% where
%   r is a structure array of rays
%   l is a scalar specifying the propagation distance
% and
%   hline is a matrix of line handles (has the same shape as r).
%
% See also: PLOT_RAYS_TO_DEST.

% Copyright 2007 by Timothy E. Holy

  if (nargin < 3)
    lineoptions = {};
  end
  hline = zeros(size(r));
  for i = 1:numel(r)
    x0 = r(i).x0;
    x1 = x0 + l*r(i).e;
    hline(i) = line([x0(1) x1(1)],[x0(2) x1(2)],'Color',r(i).rgb);
  end
  