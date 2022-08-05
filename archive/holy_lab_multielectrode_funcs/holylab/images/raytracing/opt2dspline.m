function [t_out,ro] = opt2dspline(params,r_all)
% OPT2DSPLINE: define an arbitrary interface for ray tracing
% Syntax:
%   hline = opt2dspline(params)
%   [t,ro] = opt2dspline(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces some rays.
%
% Params is the parameters structure describing the interface:
%   center: a 2-vector describing the center of the surface
%   normal: a unit 2-vector describing the direction of the "optic axis"
%     through this surface (obviously not nec. symmetric around it,
%     though)
%   xg: the "grid" (perhaps irregular) of points perpendicular to the
%     optic axis at which the surface displacement will be defined
%   yg: the displacements along the optic axis at the grid points
%   fields "mat1","mat2" are the material names (e.g., 'bk7') on opposite
%     sides of the surface.  mat1 is the input side if the dot product
%     (normal,r.e) is positive.
% rays is a structure array of input rays (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro contains the output rays after intersection.
%
% This is similar to OPT2DDMSPLINE, except this function does not assume
% that the deformations are small. It also defines a transmissive interface
% rather than a reflective surface.
%
% See also: OPT2DDMSPLINE.

% Copyright 2006 by Timothy E. Holy

if (nargin == 1)
  % Plot surface on screen
  y = linspace(params.xg(1),params.xg(end),1000);
  x = spline(params.xg,params.yg,y);
  R = o2ds_rotmtrx(params.normal);
  X = R*[x; y] + repmat(params.center(:),1,length(x));
  line(X(1,:),X(2,:),'Color','k');
else
  % Trace the rays
  n_rays = length(r_all);
  ro = r_all;
  R = o2ds_rotmtrx(params.normal);
  x0 = [r_all.x0];
  e = [r_all.e];
  x0 = x0 - repmat(params.center(:),1,n_rays);
  x0 = R*x0;  % rotate so optic axis is pointing to +x (to the right)
  e = R*e;
  yrange = [min(params.xg), max(params.xg)];
  surffunc = ppcreate(params.xg,params.yg,'spline');
  surfslope = ppcreate(surffunc,'diff');
  xstrike = nan(2,n_rays);
  strikenormal = nan(2,n_rays);
  for rayIndex = 1:n_rays
    % Check to see if ray strikes in region defined by the grid
    t0 = -x0(1,rayIndex)/e(1,rayIndex);
    y0 = x0(2,rayIndex) + t0*e(2,rayIndex);
    if (false && (y0 < yrange(1) || y0 > yrange(2)))  % fixme: the false part is a hack
      ro(rayIndex).valid = false;
      t_out(rayIndex) = t0;
      xstrike(:,rayIndex) = [0;y0];
    else
      intersect_fcn = @(t) o2ds_intersect(t,x0(:,rayIndex),...
        e(:,rayIndex),surffunc);
      t = fsolve(intersect_fcn,t0,optimset('Display','none'));
      t_out(rayIndex) = t;
      xstrike(:,rayIndex) = x0(:,rayIndex) + t*e(:,rayIndex);
      slope = surfslope(xstrike(2,rayIndex));
      strikenormal(:,rayIndex) = [1;-slope]/sqrt(1+slope^2);
    end
  end
  %xend = xstrike+strikenormal;
  %line([xstrike(1,:); xend(1,:)],[xstrike(2,:); xend(2,:)],'Color','m');
  % Now convert back to original coordinate system
  xstrike = R'*xstrike + repmat(params.center(:),1,n_rays);
  strikenormal = R'*strikenormal;
  for rayIndex = 1:n_rays
    ro(rayIndex).x0 = xstrike(:,rayIndex);
  end
  validFlag = [ro.valid];
  % Refract the rays at the surface
  ro(validFlag) = transmitted_ray(ro(validFlag),strikenormal,...
    params.mat1,params.mat2);
end  
  
function R = o2ds_rotmtrx(normal)
  R = [normal(1) normal(2); -normal(2) normal(1)]; % Rotation matrix
  % or should we take the transpose of this?
  %R = R * [0 1; -1 0]; % To compensate for the fact that y = f(x), but
                       % theta=0 corresponds to propagation in the x
                       % direction

function dy = o2ds_intersect(t,dx,raye,func)
  dy = dx(1)+t*raye(1) - func(dx(2)+t*raye(2));
