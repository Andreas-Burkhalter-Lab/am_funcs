function [r,indices] = create_inffoc_rays_for_aberrations(Ra,theta,lambda)
% CREATE_INFFOC_RAYS_FOR_ABERRATIONS: create rays for testing infinity-focused optics
% Syntax:
%   [r,indices] =  create_inffoc_rays_for_aberrations(Ra,theta,lambda)
% where
%   Ra is the radius of the back aperture (approx f*NA, where f = focal
%     length; more precisely, it's f*tan(phi), where phi is the
%     acceptance cone of the objective from an axial object in the focal
%     plane)
%   theta is the angle of the infinity-focused rays that corresponds to
%     the off-axis point (given by y/f, where y is the distance of the
%     off-axis point from the optic axis, i.e., the field of view, and f
%     is the focal length)
%   lambda is a vector of wavelengths to test the system with. The first
%     one will be the central "design wavelength" at which all
%     aberrations are measured.
% and
%   r is a structure array of rays, emanating from the z=0 plane (the
%     aperture is in this plane, too)
%   indices is a structure which gives a "name" to the various rays in r,
%     of use in functions like CALCULATE_ABERRATIONS.
%
% See also: CALCULATE_ABERRATIONS.
  
% Copyright 2007 by Timothy E. Holy
  
  r = ray;
  paraxial = 0.01*Ra;
  r.x0 = [0;0;0];
  r.e = [0;0;1];
  r.w = lambda(1);
  indices.axial = 1;
  % Create a meridional "paraxial" ray
  indices.paraxial = 2;
  i = indices.paraxial;
  r(i) = r(indices.axial);
  r(i).x0 = [0;paraxial; 0];
  % Create a meridional marginal ray
  indices.marginal = 3;
  i = indices.marginal;
  r(i) = r(indices.axial);
  r(i).x0 = [0; Ra; 0];
  % Create a meridional "zonal" ray
  indices.zonal = 4;
  i = indices.zonal;
  r(i) = r(indices.axial);
  r(i).x0 = [0; 0.7 * Ra; 0];
  % Create meridional rays corresponding to an off-axis point (so they
  % have non-zero slope)
  indices.principal = 5;
  i = indices.principal;
  r(i) = r(indices.axial);
  r(i).e = [0; sin(theta); cos(theta)];
  indices.coma = [6 7];
  i = indices.coma;
  r(i) = repmat(r(indices.principal),1,length(i));
  r(i(1)).x0 = [0; Ra; 0];
  r(i(2)).x0 = [0; -Ra; 0];
  len = length(r);
  % Rays for field curvature
  indices.tangential = len+1;
  i = indices.tangential;
  r(i) = r(indices.principal);
  r(i).x0 = [0; paraxial; 0];
  indices.sagittal = len+2;
  i = indices.sagittal;
  r(i) = r(indices.principal);
  r(i).x0 = [paraxial; 0; 0];
  len = length(r);
  % Rays for chromatic aberration
  for i = 2:length(lambda)
    r(len-1+i) = r(indices.paraxial);
    r(len-1+i).w = lambda(i);
  end
  indices.chromatic_axial = len+1:len-1+length(lambda);
  len = length(r);
  % Rays for spherochromatism
  for i = 2:length(lambda)
    r(len-1+i) = r(indices.marginal);
    r(len-1+i).w = lambda(i);
  end
  indices.spherochromatism = len+1:len-1+length(lambda);
  len = length(r);
  % Go through and provide appropriate RGB color for each wavelength
  for i = 1:length(r)
    r(i).rgb = spectrumRGB(r(i).w);
  end
  