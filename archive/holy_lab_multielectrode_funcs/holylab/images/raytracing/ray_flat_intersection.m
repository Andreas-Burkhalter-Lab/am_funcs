function [t,strike_pos,e_dot_n] = ray_flat_intersection(params,r)
% ray_flat_intersection: compute intersection of rays with flat surface
% Syntax:
%   [t,strike_pos] = ray_flat_intersection(params,rays)
%   [t,strike_pos,e_dot_n] = ray_flat_intersection(params,rays)
%
% Params is the parameters structure describing the flat surface:
%   center is the position of the center of the flat surface;
%   normal is a unit vector perpendicular to this surface;
% rays is a structure array containing info about the input rays (see RAY)
% and
%   t is the distance along the ray to the point of intersection
%   strike_pos is the position of their intersection.
%   e_dot_n is the dot product of the ray direction vector with the
%     normal (i.e., the cosine of the angle between them)
%
% See also: RAY.

  e = [r.e];
  x0 = [r.x0];
  [n_dims n_rays] = size(e);
  repcenter = repmat(params.center(:),1,n_rays);
  dx = repcenter-x0;
  N = diag(params.normal);
  e_dot_n = sum(N*e);
  t = sum(N*dx)./e_dot_n;
  strike_pos = x0 + repmat(t,n_dims,1).*e;
