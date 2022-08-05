function rb = optical_collimated_bundle(x0,e0,wavelength,outerdiameter,innerdiameter,theta,ngrid)
% optical_collimated_bundle: build a collimated beam of rays
%
% Syntax:
%   rb = optical_collimated_bundle(x0,e0,wavelength,outerdiameter,innerdiameter,theta,ngrid)
% where
%   x0 is a 3-vector giving the coordinates of the beam center
%   e0 is a 3-vector giving the direction of propagation
%   wavelength is a scalar (in nanometers)
%   outerdiameter is the total diameter of the beam
%   innerdiameter, if present and not empty, is the diameter that is within
%     1/e^2 of the peak power (what would normally be called the
%     "diameter" for a Gaussian beam). If not supplied, the beam has
%     uniform intensity.
%   theta, if supplied, causes this to produce meridional rays in the plane
%     tilted at an angle theta from the first coordinate.
%   ngrid (default 21) is the number of rays along each axis of the beam
% and
%   rb is an object of class opticalRayBundle.
%
% See also: opticalRayBundle, optical_point_bundle.

% Copyright 2010 by Timothy E. Holy

  %% Parse inputs
  if (numel(e0) ~= 3)
    error('For a collimated beam, e0 must be a single vector');
  end
  if (numel(x0) ~= 3)
    error('The origin of the collimated beam must be a single point');
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
    innerdiameter = [];
  end
  
  %% Make a bundle propagating down the z-axis starting at 0
  % Normalize e0
  e0 = e0(:) / sqrt(sum(e0.^2));
  x = linspace(-0.5,0.5,ngrid)*outerdiameter;
  if (isempty(theta))
    % A 2-dimensional grid of rays
    [X,Y] = ndgrid(x,x);
    R2 = X.^2 + Y.^2;
    keepFlag = R2 <= outerdiameter^2/4;
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
  Z = zeros(size(X));
  n_rays = size(X,1);
  
  %% Rotate so that the source plane is perpendicular to e0
  q = multiply_quaternion([0 0 0 1]',[0; e0(:)]);
  p = rotate_quaternion([X Y Z]',q);
  
  %% Displace so center is at desired location
  p = bsxfun(@plus,p,x0(:));
  
  %% Set the intensity
  if (~isempty(innerdiameter))
    I = exp(-R2(keepFlag)/(2*(innerdiameter/2)^2));
  else
    I = ones(1,n_rays);
  end
  
  %% Construct the bundle
  rb = opticalRayBundle(p,repmat(e0,1,n_rays),repmat(wavelength,1,n_rays),I);
  
  