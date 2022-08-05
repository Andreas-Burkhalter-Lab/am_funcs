function f = make_gaussian_filter(imsz,pix2phys,resolution_phys)
% make_gaussian_filter: create a spatially-transformed Gaussian
%
% This function creates a general zero-centered Gaussian suitable for use
% with Fourier filtering. The aspect ratio of the Gaussian can be
% controlled arbitrarily. The syntax is designed for supporting tilted
% imaging, like HS-OCPI configurations.
%
% Syntax:
%   f = make_gaussian_filter(imsz,pix2phys,resolution_phys)
% where
%   imsz is a vector containing the dimensions of the array that you want
%     to filter.
%   pix2phys is a matrix that expresses the displacement, in physical units
%     (e.g., microns), that arise when you move one pixel (array grid
%     spacing) along each coordinate. Making this a generic matrix allows
%     you to use non-orthogonal array coordinates.
%   resolution_phys is a vector that contains the width of the gaussian
%     along each "physical" axis.
%
% On output,
%   f is an array of size imsz, with the Gaussian parametrized in array
%     (pixel) coordinates.
%
% Example:
% Suppose pixels in y are spaced 0.3 microns apart, and the z axis
%   corresponds to a displacement dz at an angle of 30 degrees relative
%   to y. Then pix2phys is
%      theta = pi/6; dy = 0.3; pix2phys = [dy 0; -dz*cos(theta) -dz*sin(theta)];
% Suppose the sigma of resolution along y is dy, and along z it is 2.5.
% Then create f with
%    f = make_gaussian_filter(size(yz),pix2phys,[dy 2.5]);
% Display f to make sure it is close to the PSF.  
%
% You can use this filter as follows:
%   imf = ifftn(fftn(im).*fftn(f));
% You can deconvolve in the following way:
%   ffft = fftn(f);
%   finvfft = conj(ffft)./(gamma^2 + conj(ffft).*ffft);  % Wiener filter
%   imdeconv = ifftn(fftn(im).*finvfft);
% gamma is a constant that helps control how seriously to take frequencies
% with very little power. The noisier your images are, the larger gamma
% should be.  Alternatively, try
%   imdeconv = nonnegdeconvmd(im,f);
% This is slower but imposes a nonnegativity constraint, which can improve
% the results.
%
% See also: nonnegdeconvmd.

% Copyright 2011 by Timothy E. Holy

  n_dims = length(imsz);
  if length(resolution_phys) ~= n_dims
    error('The resolution must be a vector of length n_dims');
  end
  % Create the pixel coordinate system
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    tmp = 0:imsz(dimIndex)-1;
    wrapFlag = tmp > imsz(dimIndex)/2;
    tmp(wrapFlag) = tmp(wrapFlag)-imsz(dimIndex);
    x{dimIndex} = tmp;
  end
  X = cell(1,n_dims);
  [X{:}] = ndgrid(x{:});
  % Turn coordinates into a matrix
  for dimIndex = 1:n_dims
    X{dimIndex} = X{dimIndex}(:);
  end
  Xm = cat(2,X{:})';
  % Create the pixelwise resolution matrix
  A = pix2phys * diag(1./resolution_phys(:).^2) * pix2phys';
  % Calculate the filter 
  f = exp(-sum(Xm .* (A*Xm),1)/2) / (2*pi)^(n_dims/2) * sqrt(det(A));
  % Shape the filter appropriately
  f = reshape(f,imsz)/sum(f);
  