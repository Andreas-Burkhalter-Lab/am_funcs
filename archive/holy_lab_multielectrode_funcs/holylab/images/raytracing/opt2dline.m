function [t_all,ro] = opt2dline(params,r_all)
% OPT2DLINE: define a planar "surface" (line) for ray tracing
% Syntax:
%   hline = opt2dline(params)
%   [t,ro] = opt2dline(params,rays)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%
% Params is the parameters structure describing the planar surface:
%   center is the aperture center (a 2-vector)
%   normal is a unit vector perpendicular to the line
%   apr is the aperture radius
%   fields "mat1","mat2" are the material names (e.g., 'bk7') opposite
%     sides of the surface.  mat1 is the input side if the dot product
%     (normal,r.e) is positive.
% rays is the input structure array of rays (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro contains the output rays after intersection.
%
% See also: RAY, OPT2DCIRCLE.

  if (nargin == 1)
    % Plot surface on screen
    parallel = [params.normal(2);-params.normal(1)];
    xl = params.center + params.apr*parallel;
    xr = params.center - params.apr*parallel;
    t = line([xl(1) xr(1)],[xl(2) xr(2)],'Color','k');
    if isfield(params,'tag')
      set(t,'Tag',params.tag);
    end
  else
    t_all = nan(1,length(r_all));
    ro = r_all;
    validFlag = [r_all.valid];
    validIndex = find(validFlag);
    [t,strike_pos] = ray_flat_intersection(params,r_all(validFlag));
    n_rays = length(t);
    strike_dx = strike_pos - repmat(params.center(:),1,n_rays);
    % Did it strike within the aperture?
    goodFlag = sum(strike_dx.^2) <= params.apr^2;
    if any(goodFlag)
      % These are the rays that pass through the aperture
      goodIndex = validIndex(goodFlag);
      t_all(goodIndex) = t(goodFlag);
      for i = 1:length(goodIndex)
        ro(goodIndex(i)).x0 = strike_pos(:,goodIndex(i));
      end
      ro(goodIndex) = transmitted_ray(ro(goodIndex),params.normal,params.mat1,params.mat2);
    end
    if ~all(goodFlag)
      badIndex = validIndex(~goodFlag);
      [ro(badIndex).valid] = deal(false);
    end
  end
