function [xo,eo,s] = trace_rays_sineobj(x,e,f)
% TRACE_RAYS_SINEOBJ: ray-tracing through an objective lens that satisfies the sine condition.
% Syntax:
%   [xo,eo,s] = trace_rays_sineobj(x,e,f)
% where
%   x is a 3-by-n_rays matrix containing a point on each ray. The optic
%     axis is along the z-axis, and the focal plane is defined by
%     x(3,i) = 0. The sign convention is such that x(3,i) < 0 corresponds
%     to a point that is farther from the objective than the focal plane.
%   e is a 3-by-n_rays matrix, where each column is a unit vector giving
%     the propagation direction
%   f is a 2-vector, [ffront fback], of focal lengths for the lens. (See
%     GENERATE_RAYS for help in specifying this quantity.)
% and
%   xo is the set of points in the back pupil plane for each ray (a
%     3-by-n_rays matrix)
%   eo is the matrix containing the output propagation directions
%   s is a vector of optical path lengths traversed by each ray to get to
%     the output point xo.
%
% See also: TRACE_RAYS_SINEOBJTUBE, GENERATE_RAYS.

% Copyright 2008 by Timothy E. Holy

  % First trace (or back-trace) to focal plane
  dfp = -x(3,:) ./ e(3,:);
  x = x + repmat(dfp,3,1) .* e;
  
  % Trace through the objective
  eotilde = -x(1:2,:)/f(2);
  xotilde = x(1:2,:) - f(1)*e(1:2,:);
  xo = [xotilde; zeros(1,size(x,2))];
  eo = [eotilde; sqrt(1-sum(eotilde.^2,1))];
  if (nargout > 2)
    s = dfp + sum(xo.*eo,1) - f(1);
  end
end
