function [t_all,ro] = opt2dsinecond(params,r_all)
% OPT2DSINECOND: trace rays as required by the sine condition
% Syntax:
%   hline = opt2dsinecond(params)
%   [t,ro] = opt2dsinecond(params,rays)
% The first syntax plots the "lens element" on the current axis. The second
% traces some rays.
%
% This is an alternative to the standard projective-transformation formalism
% that describes perfect imaging (an unachievable goal using real optics).
% This instead describes optics that satisfy the sine condition, which is as
% close to perfect as real optics can get.  This function allows one to
% trace rays for a lens satisfying the sine condition without having to
% explicitly design the real-world elements that create such a lens.
%
% Params is the parameters structure describing the optics of the
% transformation:
%   center: a 2-vector specifying the intersection of the optic axis with
%     the vertex plane.
%   normal is the 2-vector defining the optic axis.
%   f is a 2-vector of focal lengths, [focallength1 focallength2], where
%     focallength1 is the separation from the first focal point to the
%     vertex.  These follow Born & Wolf sign convention, so that f(1)
%     measures the distance between focalpoints(:,i) and the vertex plane
%     along the (signed) optic axis.  In particular, a convergent lens has
%     focallength1 > 0, and focallength2 < 0.  Furthermore,
%           focallength1/focallength2 = -n1/n2,
%     where ni is the index of refraction in region i.
%     You can also supply a scalar for this value, in which case
%     focallength2 is set to be -focallength1.
%   thetamax (optional) if present, rejects rays that make an angle greater
%     than this value with respect to the normal. I.e., this allows you
%     to set the numerical aperture for the lens.  This can be a 2-vector
%     if you want to specify different values for the two sides.
%   plotsize is the height of the line or box used to represent this
%     element when plotted on the screen (this has no functional
%     consequences, only for visual representation; default 1)
% rays is a structure array containing info about the input rays (see RAY)
% and
%   t is the distance along the ray to the point of intersection
%     (NaN if no intersection);
%   ro contains the output rays after intersection.
%
% This is based on Meijere & Velzel, "Linear ray-propagation models in
% geometrical optics." J. Opt. Soc. Am. A 4: 2162-2165, 1987.  One point
% to note is that we use the opposite sign convention for the focal
% lengths, to be compatible with OPT2DPROJECTIVE.
%
% See also: OPT2DPROJECTIVE, RAY.

% Copyright 2007 by Timothy E. Holy

  f = (params.f(:))';  % make it a row vector
  if (length(f) == 1)
    f(2) = -f(1);  % assume a lens, rather than a mirror, unless specified
  end
  normal = params.normal(:);
  parallel = [-normal(2);normal(1)]; % parallel to lens surface
  if (nargin == 1)
    % Plot element on screen
    if ~isfield(params,'plotsize')
        params.plotsize = 1;
    end
    parallel_draw = params.plotsize*parallel;
    center = params.center;
    yp = center + parallel_draw;
    ym = center - parallel_draw;
    t = line([yp(1) ym(1)],[yp(2) ym(2)],'Color','k','LineWidth',2);
  else
    n_rays = length(r_all);
    ro = r_all;
    t_all = nan(1,n_rays);
    % Determine whether the rays are propagating along, or opposite, the
    % optic axis
    e = [r_all.e];
    N = diag(normal);
    e_dot_n = sum(N*e);
    entry_index = 2 - (e_dot_n > 0);% 1 if entering side 1, 2 if side 2
    entry_sign = 2*(e_dot_n>0) - 1; % +1 if entering side 1, -1 if side 2
    x0 = [r_all.x0];
    if isfield(params,'thetamax')
      % Test to see whether some rays fail to make it through the element
      if (length(params.thetamax) < 2)
        params.thetamax(2) = params.thetamax(1);
      end
      rejectedFlag = entry_sign.*e_dot_n < cos(params.thetamax(entry_index));
      if any(rejectedFlag)
        [ro(rejectedFlag).valid] = deal(false);
        keptFlag = ~rejectedFlag;
        e = e(:,keptFlag);
        x0 = x0(:,keptFlag);
        e_dot_n = e_dot_n(keptFlag);
        entry_index = entry_index(keptFlag);
        entry_sign = entry_sign(keptFlag);
        n_rays = length(e_dot_n);
        keptIndex = find(keptFlag);
      else
        keptIndex = 1:n_rays;
      end
    else
      keptIndex = 1:n_rays;
    end
    % Now trace the ray
    f = -f;  % Switch to M&V sign convention
    f_enter = f(entry_index);
    f_exit = f(3-entry_index);
    dx = x0 - repmat(params.center(:),1,n_rays);
    z = sum(N*dx,1) - f_enter;
    t = abs(z + f_enter); % distance to incoming circular surface
    P = diag(parallel);
    y = sum(P*dx); % dot product with vector parallel to lens
    M = sum(P*e);
    Mp = (z.*M - y)./f_exit;
    yp = y - (f_enter + z).*M; %pupil coordinate
    % Compute the displacement of the pupil exit position along the optic
    % axis. Note this fixes a sign error in Eq. 8 of M&V
    zp = (yp + f_enter.*M)./Mp .* sqrt(1-Mp.^2);
    zp(Mp == 0) = -f_exit(Mp == 0);
%     dz = f(3-entry_index) .* (z+f(entry_index)) .* Mp.^2./(1+sqrt(1-Mp.^2));
%     isz = (z == 0);
%     dz(~isz) = dz(~isz) ./ z(~isz);
%     dz(isz) = 0;
    % Take these scalar coordinates and go back to a vector
    % representation of the ray
%     x0_out = repmat(params.center(:),1,n_rays) + parallel*yp + ...
%       normal*dz;
    x0_out = repmat(params.center(:),1,n_rays) + parallel*yp + ...
      normal*(zp + f_exit);
    e_out = parallel*Mp + normal*(sqrt(1-Mp.^2).*entry_sign);
    for i = 1:length(keptIndex)
      ro(keptIndex(i)).x0 = x0_out(:,i);
      ro(keptIndex(i)).e = e_out(:,i);
    end
    t_all(keptIndex) = t;
  end
  
