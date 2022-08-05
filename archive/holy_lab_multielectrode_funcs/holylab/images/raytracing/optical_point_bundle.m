function rb = optical_point_bundle(x0,e0,wavelength,outerangle,innerangle,theta,ngrid)
% optical_point_bundle: build a beam of rays emerging from a single point
%
% Syntax:
%   rb = optical_point_bundle(x0,e0,wavelength,outerangle,innerangle,theta,ngrid)
% where
%   x0 is a 3-vector giving the coordinates of the point
%   e0 is a 3-vector giving the direction of propagation of the center ray
%   wavelength is a scalar (in nanometers)
%   outerangle is the maximum angle of the cone of rays relative to the
%     center ray. Note this might be called the "half angle" rather than
%     the "full angle."  With "hard" apertures, this is closely related to
%     the NA.
%   innerangle, if present and not empty, is the cone angle at which the
%     power becomes 1/e^2 of the peak power (what would normally be called
%     the "divergence angle" for a Gaussian beam, related to the NA of the
%     Gaussian beam). If not supplied, the beam has uniform intensity.
%   theta, if supplied, causes this to produce meridional rays in the plane
%     tilted at an angle theta from the first coordinate.
%   ngrid (default 21) is the number of rays along each axis of the beam
% and
%   rb is an object of class opticalRayBundle. The optic axis is initialized to
%     be equal to e0.
%
% See also: opticalRayBundle, optical_collimated_bundle.

% Copyright 2010 by Timothy E. Holy

  %% Parse inputs
  if (numel(e0) ~= 3)
    error('For a point beam, e0 must be a single vector');
  end
  if (numel(x0) ~= 3)
    error('The origin of the point beam must be a single point');
  end
  if ~isscalar(wavelength)
    error('wavelength must be a scalar');
  end
  if (nargin < 7 || isempty(ngrid))
    ngrid = 21;
  end
  if (nargin < 6)
    theta = [];
  end
  if (nargin < 5)
    innerangle = [];
  end
  
  %% Make a point beam symmetric around the z-axis starting at 0
  % Normalize e0
  e0 = e0(:) / sqrt(sum(e0.^2));
  soa = sin(outerangle);
  x = linspace(-1,1,ngrid)*soa;
  if (isempty(theta))
    % A 2-dimensional grid of rays
    [X,Y] = ndgrid(x,x);
    R2 = X.^2 + Y.^2;
    keepFlag = R2 <= soa^2;
    X = X(keepFlag);
    Y = Y(keepFlag);
  else
    % Rays arranged in a line
    X = x(:)*cos(theta);
    Y = x(:)*sin(theta);
    % Next two are needed in case of Gaussian profile
    R2 = X.^2 + Y.^2;
    keepFlag = true(size(R2));
  end
  X = cos(X);
  Y = sin(Y);
  Z = sqrt(X.^2 + Y.^2);
  n_rays = size(X,1);
  
  %% Rotate so that the propagation direction is along e0
  q = multiply_quaternion([0 0 0 1]',[0; e0(:)]);
  e = rotate_quaternion([X Y Z]',q);
  
  %% Set the intensity
  if (~isempty(innerangle))
    I = exp(-R2(keepFlag)/(2*innerangle^2));
  else
    I = ones(1,n_rays);
  end
  
  %% Construct the bundle
  rb = opticalRayBundle(repmat(x0(:),1,n_rays),e,repmat(wavelength,1,n_rays),I);
  rb.q = [0; e0];
  
  