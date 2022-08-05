function [t,ro] = planar(params,r)
% PLANAR: define a planar surface for ray tracing
% Syntax:
%   hsurf = planar(params)
%   [t,ro] = planar(params,ray)
% The first syntax plots the surface on the current axis. The second
% traces a ray.
%
% Params is the parameters structure describing the planar surface:
%   fields "c" and "normal" define the planar surface: x.normal = c
%   fields "mat1","mat2" are the material names (e.g., 'bk7') opposite
%     sides of the surface.  mat1 is the input side if the dot product
%     (normal,r.e) is positive.
%   apc is the aperture center (a 3-vector)
%   apr is the aperture radius, or alternatively specify rectangular
%     dimensions by ape (unit vectors along axes of aperture) and apl
%     (lengths).
% r is the input ray (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, SPHERICAL, CYLINDRICAL.

  if (nargin == 1)
    % Plot on screen
    % Use mesh?
  else
    e_norm = sum(r.e.*params.normal);
    t = (params.c-sum(r.x0.*params.normal))/e_norm;
    
    % Now create the new ray
    ro = transmitted_ray(r,t,params.normal,params.mat1,params.mat2);
  end