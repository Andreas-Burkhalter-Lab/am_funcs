function [t_out,ro] = opt2daspheric(params,r_all)
% OPT2DASPHERIC: define an aspheric "surface" (line) for ray tracing
% Syntax:
%   hline = opt2daspheric(params)
%   [t,ro] = opt2daspheric(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces some rays.
%
% Params is the parameters structure describing the planar surface:
%   center: center of the aspheric
%   optax: a unit vector indicating the orientation of the optic axis
%   func: a function handle, [z,normal] = func(x,funcparams) gives the
%     surface position and normal at a particular x (x = dist from optic
%     axis). If empty, func defaults to a subfunction in this file,
%     aspher_eq.
%   funcparams: a structure containing all the parameters needed to
%     evaluate func. There must be 2 particular fields: 'k' and 'curv'.
%   xmax: the maximum absolute value of x
%   fields "mat1","mat2" are the material names (e.g., 'bk7') on opposite
%     sides of the surface.  mat1 is the input side if the dot product
%     (Rvec,r.e) is positive.
% rays is a structure array of input rays (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro contains the output rays after intersection.
%
% See also: ASPHER_EQ, RAY, OPT2DLINE, OPT2DCIRCLE.

  if (~isfield(params,'func') || isempty(params.func))
    params.func = @aspher_eq;
  end
  if (nargin == 1)
    % Plot surface on screen
    % Calculate the surface in default orientation
    x = linspace(-params.xmax,params.xmax,101);
    y = zeros(size(x));
    for i = 1:length(x)
      y(i) = params.func(x(i),params.funcparams);
    end
    R = o2da_rotmtrx(params);
    X = R * [x;y];
    % Rotate & translate surface to actual position & orientation
    t = line(X(1,:) + params.center(1),X(2,:) + params.center(2),...
             'Color','k');
  else
    n_rays = length(r_all);
    ro = r_all;
    R = o2da_rotmtrx(params);
    for rayIndex = 1:n_rays
      r = r_all(rayIndex);
      % Calculate the intersection and refracted ray
      dx = r.x0 - params.center(:);
      % Convert to vertical optic axis (we'll convert back later)
      dx = R\dx;
      raye = R\r.e;
      % Find the point of intersection
      % First get a good starting guess analytically, by keeping only the
      % elliptical terms in the aspheric lens equation.
      k = params.funcparams.k;
      curv = params.funcparams.curv;
      a = curv * (raye(1)^2 + (k+1)*raye(2)^2);
      b = curv * (raye(1)*dx(1) + (k+1)*raye(2)*dx(2)) - raye(2);
      c = curv * (dx(1)^2 + (k+1)*dx(2)^2) - 2*dx(2);
      t0 = (-b + sqrt(b^2 - a*c))/a;  % Quadratic eq. with b->2*b
      % Now refine this guess by solving the equation numerically
      %fsolveops = optimset('Display','final');
      fsolveops = optimset('Display','off');
      [t,fval,exitflag,output] = fsolve(@(ts) o2da_intersect(...
	  ts,dx,raye,params.func,params.funcparams),t0,fsolveops);
      % Now create the new ray
      xp = dx + t*raye; % Intersection point
      tout(rayIndex) = t;
      ro(rayIndex).x0 = xp;
      ro(rayIndex).e = raye;
      [y,normal(:,rayIndex)] = params.func(xp(1),params.funcparams);
    end
    ro = transmitted_ray(ro,normal,params.mat1,params.mat2);
    % Rotate and translate back
    for rayIndex = 1:n_rays
      ro(rayIndex).e = R*ro(rayIndex).e;
      ro(rayIndex).x0 = R*ro(rayIndex).x0 + params.center(:);
    end
  
  
function R = o2da_rotmtrx(params)
  R = [params.optax; -params.optax(2) params.optax(1)]; % Rotation matrix
  % or should we take the transpose of this?
  R = R * [0 1; -1 0]; % To compensate for the fact that y = f(x), but
                       % theta=0 corresponds to propagation in the x
                       % direction

function dy = o2da_intersect(t,dx,raye,aspherefunc,aspherefuncparams)
  dy = dx(2)+t*raye(2) - aspherefunc(dx(1)+t*raye(1),aspherefuncparams);
