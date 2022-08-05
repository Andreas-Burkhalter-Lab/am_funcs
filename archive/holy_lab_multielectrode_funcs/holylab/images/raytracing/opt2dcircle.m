function [t,ro] = opt2dcircle(params,r_all)
% OPT2DLINE: define a circular "surface" (line) for ray tracing
% Syntax:
%   hline = opt2dcircle(params)
%   [t,ro] = opt2dcircle(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces some rays.
%
% Params is the parameters structure describing the circular surface:
%   center: center of the circle
%   Rvec: a vector whose magnitude is the radius of the circle, and which
%     points in the direction of the center of the arc
%   theta: the angle from the center of the arc to its edge
%   fields "mat1","mat2" are the material names (e.g., 'bk7') on opposite
%     sides of the surface.  mat1 is the input side if the dot product
%     (Rvec,r.e) is positive.
% r is an input ray structure array (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, OPT2DLINE.

  if (nargin == 1)
    % Plot surface on screen
    thetac = atan2(params.Rvec(2),params.Rvec(1)); % Angle of arc center
    theta = linspace(thetac-params.theta,thetac+params.theta,100);
    R = sqrt(sum(params.Rvec .* params.Rvec));
    t = line(params.center(1)+R*cos(theta),params.center(2)+R*sin(theta),...
             'Color','k');
  else
    n_rays = length(r_all);
    ro = r_all;
    for rayIndex = 1:n_rays
      r = r_all(rayIndex);
      dx = r.x0 - params.center(:);
      dx2 = sum(dx.*dx);
      e_dx = sum(r.e.*dx);
      R2 = sum(params.Rvec.*params.Rvec);
      % Solve quadratic equation for intersection
      sqrtarg = e_dx^2 - dx2 + R2;
      if (sqrtarg < 0)
        % No intersection
        t(rayIndex) = NaN;
        ro(rayIndex).valid = false;
      end
      t(rayIndex) = -e_dx + sign(sum(r.e .* params.Rvec(:))) * sqrt(sqrtarg); %hack
      % Now create the new ray
      xp = r.x0 + r.e*t(rayIndex);
      normal = xp-params.center(:);
      normal = normal/sqrt(sum(normal.*normal));
      normal_all(:,rayIndex) = normal;
      ro(rayIndex).x0 = xp;
    end
    ro = transmitted_ray(ro,normal_all,params.mat1,params.mat2);
  end
  