function [H,rho,theta,c0] = pupil_initialize(NA,M,lambda,imsz,pixel_spacing,zstep)
% pupil_initialize: calculate parameters for Fourier optics based on camera & optical parameters
% Syntax:
%   [H,rho,theta] = pupil_initialize(NA,M,lambda,imsz,pixel_spacing)
%   [H,rho,theta,c0] = pupil_initialize(NA,M,lambda,imsz,pixel_spacing,zstep)
% where
%   NA is the numerical aperture of the objective
%   M is the magnification
%   lambda is the wavelength of light, in units of microns (or really,
%     the same units used for pixel_spacing)
%   imsz is the number of pixels along each coordinate (i.e., the size of
%     the images produced by the camera)
%   pixel_spacing should be supplied in _image_ space units, e.g., 4
%     microns for a camera with 4 micron pixels.  If your camera has
%     asymmetric pixels, supply this as a 2-vector for the pixel spacing
%     along each axis.
%   zstep (optional) is the spacing between frames in an image stack,
%     useful if you want to consider 3d imaging.  This is supplied in
%     _object_ space units (i.e., supply 1 if your images are acquired by
%     translating the objective by 1 micron).
% and
%   H is the pupil function (1 inside the pupil, 0 outside the pupil) in
%     appropriate Fourier coordinates (mapped to the same size grid as the
%     image)
%   rho is the normalized pupil radius (rho = 1 on the edge of the pupil);
%     note that rho increases beyond 1 outside the pupil
%   theta is the angular coordinate
%   c0 is the coefficient of z in setting the defocus, in units of integer
%     coordinates within the stack (i.e., the frame number)

% Copyright 2008 by Timothy E. Holy

  do3d = (nargout > 3);
  if isstruct(NA)
    % imagingConfiguration,imsz syntax
    imagingConfiguration = NA;
    imsz = M;
    NA = imagingConfiguration.NA;
    M = imagingConfiguration.M;
    lambda = imagingConfiguration.lambda;
    pixel_spacing = imagingConfiguration.pixel_spacing;
    if isfield(imagingConfiguration,'zstep')
      zstep = imagingConfiguration.zstep;
    end
  end
  
  NA_image = NA/M;
  % Set up the relationships between physical space, fourier space, and
  % the aperture function
  n_dims = 2;
  coords = cell(1,n_dims);
  for i = 1:n_dims
    coords{i} = 0:imsz(i)-1;
  end
  X = cell(1,n_dims);  % coordinates in the fourier plane
  [X{:}] = ndgrid(coords{:});
  res_factor = pixel_spacing * NA_image / lambda;
  if any(res_factor > 1)
    warning(['Pixel spacing is determining the resolution (you''re not' ...
	     ' doing diffraction-limited imaging). You will get aliasing.']);
    res_factor(res_factor > 1) = 1;
  end
  if isscalar(res_factor)
    res_factor = res_factor * [1 1];
  end
  for i = 1:n_dims
    % Set up circular boundary conditions with origin at corner
    X{i} = modwrap(X{i},imsz(i));
    % Normalize to the pupil size (so rho = 1 corresponds to pupil edge)
    X{i} = X{i} / (imsz(i) * res_factor(i));
  end
  [theta,rho] = cart2pol(X{2},X{1});
  H = double(rho <= 1);
%   h = zeros(size(rho),'single');
%   h(1,1) = 1;
%   H = fft2(h);
%   H(rho > 1) = 0;

  % Determine the coefficient that controls defocus
  % From Eqs. 8-9 of Born & Wolf, Chapter 8.8 ("The Three Dimensional Light Distribution Near Focus")  
  % (note NA = a/f on the object side)
  if do3d
    c0 = (2*pi/lambda) * NA^2 * zstep;
  end
end

function xo = modwrap(xi,modval)
% xo = xi % modval, except return in range [-modval/2,modval/2]
  xo = mod(xi,modval);
  toobig = xo >= modval/2;
  xo(toobig) = xo(toobig) - modval;
end
