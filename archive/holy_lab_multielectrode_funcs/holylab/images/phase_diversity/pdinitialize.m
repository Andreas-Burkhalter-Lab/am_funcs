function [H,pdk,rho,theta] = pdinitialize(NAred,M,lambda,imsz,pixel_spacing,zk)
% PDINITIALIZE: initialize pupil objects for phase diversity calculations
%
% This is a thin wrapper of pupil_initialize, providing slightly different
% outputs in a form useful for phase diversity calculations.
%
% Syntax:
%   [H,pdk,rho,theta] = pdinitialize(NAred,M,lambda,imsz,pixel_spacing,zk)
% where
%   NAred is reduced numerical aperture of the objective, NA/n, where n is
%     the refractive index of the immersion medium.
%   M is magnification
%   lambda is wavelength of light (typically supplied in microns, but
%     whatever units you choose for lambda need to be maintained for all
%     other measurements with units of length)
%   imsz is the number of pixels along each coordinate
%   pixel_spacing should be supplied in physical (image-space) units,
%     e.g., 4 microns for a camera with 4 micron pixels.
%   zk is the actual amount of focal displacement, measured in image-space
%     units (i.e., the amount that the camera is shifted by). This can be a
%     vector of displacements, e.g., [0 1000 2000] to denote 3 diversity
%     images: (nominal) focus, defocused by 1mm, and defocused by 2mm
%     (assuming lambda is supplied in microns).
% and
%   H is the aperture function (0 or 1)
%   pdk is a [imx imy n_images] array containing the "known" aberration for
%     each diversity image
%   rho, theta are polar coordinates in the pupil
%
% NOTE: the final input is in different format from the zstep input of
% pupil_initialize.
%
% See also: pupil_initialize.

% Copyright 2008 by Timothy E. Holy

  [H,rho,theta] = pupil_initialize(NAred,M,lambda,imsz,pixel_spacing);

  % Compute the phase diversity defocus, assuming defocus by a particular z
  % From Eqs. 8-9 of Born & Wolf, Chapter 8.8 ("The Three Dimensional Light Distribution Near Focus")  
  % (note NAred/M = a/f)
  coef = (2*pi/lambda)*(NAred/M)^2/2;
  rho2 = rho.^2;
  pdk = zeros([size(rho) length(zk)],'single');
  n_dims = length(imsz);
  colons = repmat({':'},1,n_dims);
  for i = 1:length(zk)
    pdk(colons{:},i) = rho2 * (coef * zk(i));
  end
end

