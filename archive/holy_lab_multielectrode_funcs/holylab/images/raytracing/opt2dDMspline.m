function [t_all,ro] = opt2dDMspline(params,r_all)
% OPT2DDMFLAT: a deformable mirror for ray tracing
% Syntax:
%   hline = opt2dDMspline(params)
%   [t,ro] = opt2dDMspline(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%
% This function approximates the DM as a spline surface.  For efficiency,
% this makes the approximation that the "horizontal" position (along the
% line parallel to the DM surface) of the ray strike can be calculated
% without considering the mirror's deformation. This approximation is
% most appropriate for small deformations (which should always be true) and/or
% relatively small strike angles.
%
% Params is the parameters structure describing the deformable surface:
%   center: a 2-vector of the center position
%   theta: the angle with respect to the positive-x axis of the average
%     tilt
%   length: the size of the mirror
%   DMparams: a vector of positions, describing the distance of the
%     surface perpendicular to the "plane" of the mirror. These positions
%     are evenly distributed along the length of the mirror, and the
%     position at an arbitrary point is spline-interpolated at the point
%     of ray intersection.
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
  N = diag(normal);
  c = sum(normal .* params.center(:));
  lgrid = linspace(-params.length/2,params.length/2,...
		   length(params.DMparams));
  % Calculate the ray intersections with the mirror, at first ignoring
  % the deformations
  ro = r_all;
  validFlag = [r_all.valid];
  e = [r_all(validFlag).e];
  x0 = [r_all(validFlag).x0];
  validIndex = find(validFlag);
  n_rays = size(e,2);
  e_norm = sum(N*e);
  t = (c-sum(N*x0))./e_norm;
  % Later, this will be first-order corrected to include deformations
  strike_pos0 = x0 + repmat(t,2,1).*e;
  repcenter = repmat(params.center(:),1,n_rays);
  strike_dx = strike_pos0 - repcenter;
  P = diag([cos(params.theta),sin(params.theta)]); % parallel "vector"
  strike_l = sum(P*strike_dx);  % signed distance from mirror center
  goodFlag = abs(strike_l) < params.length/2; % did it go over the edge?
  if ~all(goodFlag)
    % These are the rays that didn't hit the mirror
    invalidIndex = validIndex(~goodFlag);
    [ro(invalidIndex).valid] = deal(false);
    t_all(invalidIndex) = NaN;
  end
  if any(goodFlag)
    % For the rays that hit the mirror, calculate the mirror position and
    % slope from the deformable component
    validIndex = validIndex(goodFlag);
    mp = ppcreate(lgrid,params.DMparams,'natural'); % spline approximation
    mpd = ppcreate(mp,'diff');
    DMdist = mp(strike_l(goodFlag));
    DMslope = mpd(strike_l(goodFlag));
    % Now approximate the mirror as being piecewise linear at each
    % 0th-order strike location, and re-calculate the strike position
    % The DMslope is like tan(deltatheta), but since it will be tiny we
    % can approximate tan(deltatheta) = deltatheta
    theta = params.theta+DMslope;
    normal = [-sin(theta); cos(theta)];
    cnew = c + DMdist;
    e = [r_all(validIndex).e];
    x0 = [r_all(validIndex).x0];
    e_norm = sum(normal .* e);
    t_all(validIndex) = (c - sum(normal .* x0))./e_norm;
    strike_pos = x0 + repmat(t(validIndex),2,1) .* e;
    sth = sin(2*theta);
    cth = cos(2*theta);
    n_rays = size(strike_pos,2);
    for i = 1:n_rays
      ro(validIndex(i)).x0 = strike_pos(:,i);
      ro(validIndex(i)).e = [ cth(i) sth(i);
           sth(i) -cth(i) ] * e(:,i);
    end
  end
