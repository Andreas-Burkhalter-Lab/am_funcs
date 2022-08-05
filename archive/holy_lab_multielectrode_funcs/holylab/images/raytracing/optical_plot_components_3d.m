function hsurfo = optical_plot_components_3d(component)
% optical_plot_components_3d: draw component surfaces in three dimensions
%
% Syntax:
%   hsurf = optical_plot_components_3d(component)
% where
%   component is a cell array of materials & surfaces
% and
%   hsurf is a vector of surface handles, one per optical surface.

% Copyright 2010 by Timothy E. Holy

  surfaceFlag = cellfun(@(c) isa(c,'opticalAperture'),component);
  hsurf = zeros(1,sum(surfaceFlag));
  surfaceIndex = find(surfaceFlag);
  hold on
  for i = 1:length(surfaceIndex)
    thisIndex = surfaceIndex(i);
    [X,Y,Z] = component{thisIndex}.surf3d;
    hsurf(i) = surf(X,Y,Z);
  end
  hold off
  if (nargout > 0)
    hsurfo = hsurf;
  end
