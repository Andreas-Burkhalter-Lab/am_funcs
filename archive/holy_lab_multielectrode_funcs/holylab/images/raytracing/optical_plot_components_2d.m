function hline = optical_plot_components_2d(component,ax)
% optical_plot_components_2d: draw two-dimensional projections of components
%
% Syntax:
%   hline = optical_plot_components_2d(component,ax)
% where
%   component is a cell array of materials & surfaces
%   ax is an axis perpendicular to the optic axis in the plane of the
%     projection
% and
%   hline is a vector of line handles, one per surface.

% Copyright 2010 by Timothy E. Holy

  surfaceFlag = cellfun(@(c) isa(c,'opticalAperture'),component);
  hline = zeros(1,sum(surfaceFlag));
  surfaceIndex = find(surfaceFlag);
  for i = 1:surfaceIndex
    thisIndex = surfaceIndex(i);
    [y,z] = component{thisIndex}.surf2d(ax);
    hline(i) = line(z,y,'Color','k');
  end
