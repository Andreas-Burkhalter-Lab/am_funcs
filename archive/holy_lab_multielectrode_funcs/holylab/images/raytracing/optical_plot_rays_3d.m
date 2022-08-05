function hlineo = optical_plot_rays_3d(rb,col)
% optical_plot_rays_3d: plot traced rays
%
% Syntax:
%   hline = optical_plot_rays_3d(rb)
%   hline = optical_plot_rays_3d(rb,col)
% where
%   rb is a vector of opticalRayBundles, at successive stages of raytracing
%   col is an RGB triple (if absence, the wavelength+intensity is converted
%     into a color individually for each ray)
% and
%   hline is a vector of line handles, one per ray.
%
% See also: opticalRayBundle.

% Copyright 2010 by Timothy E. Holy

  p = cat(3,rb.p_all);
  p = permute(p,[3 2 1]);  % convert to stage, ray#, coordinate ordering
  X = p(:,:,1);
  Y = p(:,:,2);
  Z = p(:,:,3);
  n_rays = size(X,2);
  if (nargin < 2)
    % Convert the wavelength into an RGB value
    w = rb(1).wavelength;
    [uw,~,colIndex] = unique(w);
    ucol = squeeze(spectrumRGB(uw));
    col = ucol(colIndex,:);
    I = rb(1).intensity;
    if ~all(I == 1)
      col = lightencolors(col,1-I');
    end
  end
  % Draw the lines
  hline = line(X,Y,Z);
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
