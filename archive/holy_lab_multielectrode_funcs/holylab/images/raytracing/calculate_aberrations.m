function [aberrations,zf] = calculate_aberrations(r,indices)
% CALCULATE_ABERRATIONS: calculate the classical aberrations of a set of rays
% Syntax:
%   [aberrations,zf] = calculate_aberrations(r,indices)
% where
%   r is a structure array of "output" rays, likely obtained from raytracing
%     (e.g., from OPT3D_RAYTRACE_SPHERICAL_ONAXIS)
%   indices is the set of labels returned by
%     CREATE_INFFOC_RAYS_FOR_ABERRATIONS
% and
%   aberrations is a structure whose fields list the values of the major
%     aberrations
%   zf is the z-position of the focal plane
%
% Adapted from Chapter 10 of Warren J. Smith's "Modern Optical Engineering"
%
% See also: OPT3D_RAYTRACE_SPHERICAL_ONAXIS, CREATE_INFFOC_RAYS_FOR_ABERRATIONS.
% Copyright 2007 by Timothy E. Holy
  
  if (size(r,1) > 1 && size(r,2) > 1)
    error('Only supply the final rays');
  end
  % Find the position of the focal plane
  zf = axial_intercept(r(indices.paraxial));
  % Find the axial intercept of the marginal & zonal rays
  aberrations.spherical_longitudinal = ...
      axial_intercept(r(indices.marginal)) - zf; 
  aberrations.zonal_longitudinal = ...
      axial_intercept(r(indices.zonal)) - zf;
  rc = r(indices.marginal);
  aberrations.spherical_transverse = ...
      - aberrations.spherical_longitudinal * rc.e(2)/rc.e(3);
  rc = r(indices.zonal);
  aberrations.zonal_transverse = ...
      - aberrations.zonal_longitudinal * rc.e(2)/rc.e(3);
  % Calculate the coma (approximate)
  n_coma = length(indices.coma);
  h = zeros(1,n_coma);
  for i = 1:n_coma
    tmp = focalplane_position(r(indices.coma(i)),zf);
    h(i) = tmp(2);
  end
  principal_position = focalplane_position(r(indices.principal),zf);
  aberrations.coma = mean(h) - principal_position(2);
  % Calculate the field curvature
  rc = r(indices.principal);
  tanP = rc.e(2)/rc.e(3);
  tmp = focalplane_position(r(indices.tangential),zf);
  rc = r(indices.tangential);
  tanT = rc.e(2)/rc.e(3);
  aberrations.curv_tangential = (principal_position(2) - tmp(2))/(tanP-tanT);
  tmp = focalplane_position(r(indices.sagittal),zf);
  rc = r(indices.sagittal);
  tanS = rc.e(2)/rc.e(3);
  aberrations.curv_sagittal = (principal_position(2) - tmp(2))/(tanP-tanS);
  % Chromatic aberration
  for i = 1:length(indices.chromatic_axial)
    ztmp = axial_intercept(r(indices.chromatic_axial(i)));
    aberrations.chromatic_axis(i) = ztmp - zf;
  end
  % Spherochromatism
  for i = 1:length(indices.spherochromatism)
    ztmp = axial_intercept(r(indices.spherochromatism(i)));
    aberrations.spherochromatism(i) = ztmp - zf ...
	- aberrations.spherical_longitudinal;
  end
  
function z = axial_intercept(r)
  exy = r.e(1:2);
  t = -r.x0(1:2)'*exy/(exy'*exy);
  z = r.x0(3) + t*r.e(3);
  
function x = focalplane_position(r,z)
  t = (z - r.x0(3))/r.e(3);
  x = r.x0 + t*r.e;
  