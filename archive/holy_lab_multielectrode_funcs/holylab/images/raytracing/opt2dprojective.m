function [t_all,ro] = opt2dprojective(params,r_all)
% OPT2DPROJECTIVE: paraxial lens approximation for all rays
% Syntax:
%   hline = opt2dprojective(params)
%   [t,ro] = opt2dprojective(params,rays)
% The first syntax plots the "lens element" on the current axis. The second
% traces some rays.
%
% Params is the parameters structure describing the projective
% transformation:
%   focalpoints: a two-by-two matrix, where the columns contain the
%     positions of the focal points.
%   f is a 2-vector of focal lengths, [focallength1 focallength2],
%     where focallength1 is the focal length corresponding to
%     focalpoints(:,1) and focallength2 is the same for focalpoints(:,2).
%     These follow Born & Wolf sign convention, so that f(1) measures the
%     distance between focalpoints(:,i) and the ith unit plane along the
%     (signed) optic axis.  In particular, a convergent lens has
%     focallength1 > 0, and focallength2 < 0.  Furthermore, 
%           focallength1/focallength2 = -n1/n2,
%     where ni is the index of refraction in region i.
%     You can also supply a scalar for this value, in which case
%     focallength2 is set to be -focallength1.
%   normal is the 2-vector defining the optic axis. This is parallel to
%     focalpoints(:,2) - focalpoints(:,1). It is used only for deciding the
%     direction that a rays enters from (i.e., if the dot product of
%     "normal" and the ray's "e" is positive, then it's assumed to enter
%     from region 1; if it's negative, then it enters from region 2).
%   thetamax (optional) if present, rejects rays that make an angle greater
%     than this value with respect to the normal. This can be a 2-vector
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
% See also: OPT2DSINECOND, OPT2DCIRCLE, RAY.

% Copyright 2006 by Timothy E. Holy

  f = (params.f(:))';  % make it a row vector
  if (length(f) == 1)
    f(2) = -f(1);  % assume a lens, rather than a mirror, unless specified
  end
  if ~isfield(params,'focalpoints')
    % Handle the "center, normal" input style
    center = params.center(:);
    params.focalpoints = center(:,[1 1]) - params.normal(:)*f;
  end
  
  fpvec = diff(params.focalpoints,1,2);
  fpsep = norm(fpvec);
  %normal = fpvec/fpsep; % unit vector along optic axis
  normal = params.normal;
  parallel = [-normal(2);normal(1)]; % parallel to lens surface
  unit_points = params.focalpoints + normal*f;
  if (nargin == 1)
    % Plot element on screen
    if ~isfield(params,'plotsize')
        params.plotsize = 1;
    end
    parallel_draw = params.plotsize*parallel;
    for i = 1:2
      center = unit_points(:,i);
      yp = center + parallel_draw;
      ym = center - parallel_draw;
      t = line([yp(1) ym(1)],[yp(2) ym(2)],'Color','k','LineWidth',2);
    end
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
    % Calculate the point of intersection of each ray
    % with the input unit "plane" (really a line, since this is 2d)
    entry_up = unit_points(:,entry_index); % the entry unit point
    x0 = [r_all.x0];
    dx = entry_up-x0;
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
        dx = dx(:,keptFlag);
        e_dot_n = e_dot_n(keptFlag);
        entry_up = entry_up(:,keptFlag);
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
    t = sum(N*dx)./e_dot_n;
    strike_pos = x0 + repmat(t,2,1).*e;
    strike_dx = strike_pos - entry_up;
    P = diag(parallel);
    h = sum(P*strike_dx); % the "height" of the strike
    % Now we apply the rules of the projective transformation to
    % calculate the new propagation direction:
    %     f' tan(gamma') + f tan(gamma) = h
    % where gamma represents the angle made with the optic axis (normal)
    tan_gamma_in = sum(P*e) ./ e_dot_n; 
    gamma_out = atan( (h - f(entry_index).*tan_gamma_in) ./ f(3-entry_index) );
    e_out = normal*(entry_sign.*cos(gamma_out)) + ...
      parallel*(entry_sign.*sin(gamma_out));
    for i = 1:n_rays
      % The ray leaves from the exit unit plane
      ro(keptIndex(i)).x0 = unit_points(:,3-entry_index(i)) + h(i)*parallel;%strike_pos(:,i);
      ro(keptIndex(i)).e = e_out(:,i);
    end
    t_all(keptIndex) = t;
  end
  
