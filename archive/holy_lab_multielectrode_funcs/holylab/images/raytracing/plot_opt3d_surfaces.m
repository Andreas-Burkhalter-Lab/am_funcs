function plot_opt3d_surfaces(surfaces,r,options)
% PLOT_OPT3D_SURFACES: display surfaces used for 3d raytracing
% Syntax:
%   plot_opt3d_surfaces(surfaces,r)
%   plot_opt3d_surfaces(surfaces,r,options)
% where
%   surfaces is of the form described in OPT3D_RAYTRACE_SPHERICAL_ONAXIS.
%   r is the output raytrace from OPT3D_RAYTRACE_SPHERICAL_ONAXIS (used to
%     determine the minimum diameter of each optical component)
%   options is a structure which may have the following fields
%     shading (default false): if true, each optical component is shaded in
%       proportion to its index of refraction (higher = darker).
%     n (required if shading is true): refractive index of each surface
%       (will be of length n_surfaces+1, specifying material on either side
%       of lens system)
%     max_darkness (default 1): from 0 = no shading to 1 = highest index
%       is black
%     coverslip (default []): if set to [height thick], the final material
%       will plotted as a rectangular surface of the given height and
%       thickness to represent a coverslip.
%
% The perspective is of the meridional plane.
%
% See also: OPT3D_RAYTRACE_SPHERICAL_ONAXIS.
  
% Copyright 2008 by Timothy E. Holy.
  
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'shading',false,'max_darkness',1,'coverslip',[]);

  z0 = 0;
  ybase = linspace(-1,1,101);
  lasty = [];
  lastz = [];
  
  objhandles = [];
  for surfIndex = 1:length(surfaces)
    % Determine the max height from the rays
    rtmp = r(:,surfIndex+1)';
    pos = [rtmp.x0];
    mh = sqrt(max(sum(pos(1:2,:).^2)));
    y = ybase*mh;
    % Calculate the z-position of the surface at each height
    c = surfaces(surfIndex).c;
    dz = c * y.^2 ./ (1 + sqrt(1 - c^2 * y.^2));
    z0 = z0 + surfaces(surfIndex).t;
    % Draw the surface
    if options.shading
      col = options.max_darkness * (options.n(surfIndex)-1)/(max(options.n)-1);
      col = (1-col) * [1 1 1];
      if (surfIndex > 1)
        combinedx = [lastz(end:-1:1) z0+dz];
        combinedy = [lasty(end:-1:1) y];
        objhandles(end+1) = patch(combinedx,combinedy,col);
      end
    else
      col = 'k';
    end
    lastz = z0+dz;
    lasty = y;
    objhandles(end+1) = line(z0+dz,y,'Color',col);
  end
  if ~isempty(options.coverslip)
    if options.shading
      col = options.max_darkness * (options.n(end)-1)/(max(options.n)-1);
      col = (1-col) * [1 1 1];
      x = lastz(end) + [0 0 1 1]*options.coverslip(2);
      y = [-1 1 1 -1]*options.coverslip(1);
      objhandles(end+1) = patch(x,y,col);
    end
  end
      
  % Put the surfaces on the bottom layer, by placing them first in the
  % drawing order
  if options.shading
    hc = get(gca,'Children');
    indx = findainb(objhandles,hc);
    hc(indx) = [];
    hc = [hc; objhandles(:)];
    set(gca,'Children',hc);
  end
  