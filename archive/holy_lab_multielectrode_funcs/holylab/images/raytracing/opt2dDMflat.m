function [t_all,ro] = opt2dDMflat(params,r_all)
% OPT2DDMFLAT: a simple deformable mirror for ray tracing
% Syntax:
%   hline = opt2dDMflat(params)
%   [t,ro] = opt2dDMflat(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%
% This function approximates the DM as a flat surface for the purpose of
% intersection with the ray; the only effect of the DM is to provide a
% change of angle to the reflected ray.  This "approximation" is
% exact for LCoS SLMs, and pretty good for mirror DMs because of
% their limited stroke.
%
% Params is the parameters structure describing the deformable surface:
%   center: a 2-vector of the center position
%   theta: the angle with respect to the positive-x axis of the average
%     tilt
%   length: the size of the mirror
%   DMparams: a vector of local tilts (deltatheta), describing the
%     warping of the surface. These are evenly distributed along the
%     length of the mirror, and the deltatheta at the point of ray
%     intersection is interpolated linearly between them.
%
% rays is a structure array of input rays (see RAY)
%
% On output,
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro contains the output rays after intersection.
%
% See also: RAY, OPT2DLINE, OPT2DCIRCLE, OPT2DMIRROR


  if (nargin == 1)
    % Plot surface on screen 
    dx = params.length*[cos(params.theta),sin(params.theta)]/2;
    line(params.center(1)+dx(1)*[-1 1],params.center(2)+dx(2)*[-1 1],...
      'Color','k','LineWidth',2);
    return
  end
  
  % Calculate alternative parametrization of mirror position
  normal = [-sin(params.theta);cos(params.theta)];
  c = sum(normal .* params.center(:));
  lgrid = linspace(-params.length/2,params.length/2,...
		   length(params.DMparams));
  % Loop over rays
  ro = r_all;
  n_rays = length(r_all);
  for rayIndex = 1:n_rays
    r = r_all(rayIndex);
    % Calculate ray intersection
    e_norm = sum(r.e .* normal);
    t = (c - sum(r.x0 .* normal))/e_norm;
    strike_pos = r.x0 + t*r.e;
    strike_dist = norm(strike_pos - params.center);
    if (strike_dist > params.length/2 || t < 0)
      % It did not hit the mirror
      ro(rayIndex).valid = false;
      t_all(rayIndex) = NaN;
      l(rayIndex) = 0;
    else
      t_all(rayIndex) = t;
      ro(rayIndex).x0 = strike_pos;
      % Set up for calculation of deltatheta at the strike position
      parallel = [normal(2);-normal(1)];
      l(rayIndex) = sum((strike_pos - params.center) .* parallel);
    end
  end
  % Now do the interpolation of the angle deltatheta
  deltatheta = interp1(lgrid,params.DMparams,l);
  % Calculate the new propagation direction
  theta = params.theta + deltatheta;
  sth = sin(2*theta);
  cth = cos(2*theta);
  for rayIndex = 1:n_rays
    ro(rayIndex).e = [ cth(rayIndex) sth(rayIndex);
           sth(rayIndex) -cth(rayIndex) ] * ro(rayIndex).e;
  end
