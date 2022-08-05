function hline = plot_rays_to_dest(r,x)
% plot_rays_to_dest: plot rays as they head towards a particular point
% Syntax:
%   hline = plot_rays_to_dest(r,x)
% where
%   r is a structure array of rays
%   x is a 2-by-1 vector indicating the point to which the ray might
%       propagate
% and
%   hline is a matrix of line handles (has the same shape as r).
%
% Note the rays to not actually have to pass through x; they are simply
% plotted for a sufficient distance to do so.
%
% See also: PLOT_RAYS_FOR_LENGTH.
  
% Copyright 2007 by Timothy E. Holy

  if (nargin < 3)
    lineoptions = {};
  end
  hline = zeros(size(r));
  for i = 1:numel(r)
    dx = x - r(i).x0;
    t = sqrt(sum(dx.^2));
    % Go 10% farther
    t = 1.1*t;
    if (sum(dx .* r(i).e) < 0)
      t = -t;  % backtrace if need be
    end
    x0 = r(i).x0;
    x1 = x0 + t*r(i).e;
    % If this is a 3d ray, assume we want a z,y plot
    if (length(x0) == 3)
      x0 = x0([3 2]);
      x1 = x1([3 2]);
    end
    hline(i) = line([x0(1) x1(1)],[x0(2) x1(2)],'Color',r(i).rgb);
  end
  