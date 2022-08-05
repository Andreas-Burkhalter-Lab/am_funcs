function hlineo = optical_plot_rays_2d(rb,theta,col)
% optical_plot_rays_2d: plot traced rays in a plane containing the optic axis
%
% Syntax:
%   hline = optical_plot_rays_2d(rb,theta)
%   hline = optical_plot_rays_2d(rb,theta,col)
% where
%   rb is a vector of optRayBundles, at successive stages of raytracing
%   theta is the angle of the plane to which to project the rays
%   col is an RGB triple (if absence, the wavelength+intensity is converted
%     into a color individually for each ray)
% and
%   hline is a vector of line handles, one per ray.
%
% See also: opticalRayBundle.

% Copyright 2010 by Timothy E. Holy

  p = cat(3,rb.p_all);
  p = permute(p,[3 2 1]);  % convert to stage, ray#, coordinate ordering
  X = p(:,:,1)*cos(theta) + p(:,:,2)*sin(theta);
  Z = p(:,:,3);
  n_rays = size(X,2);
  if (nargin < 3)
    % Convert the wavelength into an RGB value
    w = rb(1).wavelength_all;
    [uw,~,colIndex] = unique(w);
    ucol = squeeze(spectrumRGB(uw));
    col = ucol(colIndex,:);
    I = rb(1).intensity;
    if ~all(I == 1)
      col = lightencolors(col,1-I');
    end
  end
  % Draw the lines
  hline = line(Z,X);
  % Set the colors
  if (size(col,2) == 1)
    set(hline,'Color',col)
  else
    colc = mat2cell(col,ones(1,n_rays),3);
    set(hline,{'Color'},colc);
  end
  if (nargout > 0)
    hlineo = hline;
  end
