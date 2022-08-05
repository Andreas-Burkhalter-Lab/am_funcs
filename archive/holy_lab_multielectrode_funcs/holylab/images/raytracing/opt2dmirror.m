function [t_all,ro] = opt2dmirror(params,r_all)
% OPT2DMIRROR: define a flat mirror for ray tracing
% Syntax:
%   hline = opt2dmirror(params)
%   [t,ro] = opt2dmirror(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces some rays.
%
% Params is the parameters structure describing the planar surface:
%   center: middle of the 2D mirror (y axis)
%   length: center to edge length
%   theta: the angle from optical axis to mirror (in degrees)
%   mirPos: Position of mirror
% r is the input ray (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, OPT2DLINE, OPT2DCIRCLE, OPT2DDM

  if (nargin == 1)
    % Plot surface on screen   
    mirror_points = linspace(-params.length, params.length, 100);
    t = line(mirror_points * cos(params.theta) + params.mirPos, mirror_points * sin(params.theta) + params.center, 'Color','k','LineWidth',2);
  else
    ro = r_all;
    params.normal = [cos(pi/2 + params.theta) sin(pi/2 + params.theta)];     
    params.c = sum([params.mirPos params.center] .* params.normal);

    n_rays = length(r_all);
    R = [ cos(2*params.theta) sin(2*params.theta);...
           sin(2*params.theta) -cos(2*params.theta)];
    for rayIndex = 1:n_rays
      r = r_all(rayIndex);
      e_norm = sum(r.e.*params.normal);
      t = (params.c-sum(r.x0.*params.normal))/e_norm;
      ro(rayIndex).x0 = r.x0 + t*r.e;  
        
      t_all(rayIndex) = t;
      if (ro.x0(2) > (params.length*sin(params.theta) - params.center) || ...
	  ro.x0(2) < (-params.length*sin(params.theta) - params.center) || ...
	  t < 0)
        % It did not hit the mirror
        ro(rayIndex).valid = false;
        t_all(rayIndex) = NaN;
      end
      ro(rayIndex).e = R*r.e;
    end
  end
  
  % ro is output from mirror
  % t is the x co-ordinate of point of intersection of input rays and
  % mirror
  