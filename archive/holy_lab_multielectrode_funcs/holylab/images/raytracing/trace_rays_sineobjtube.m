function [xo,eo,s] = trace_rays_sineobjtube(x,e,fa,fb,separation)
% TRACE_RAYS_SINEOBJTUBE: ray-tracing through a microscope that satisfies the sine condition.
% Syntax:
%   [xo,eo,s] = trace_rays_sineobjtube(x,e,fobj,ftube,separation)
% where
%   x is a 3-by-n_rays matrix containing a point on each ray. The optic
%     axis is along the z-axis, and the focal plane is defined by
%     x(3,i) = 0. The sign convention is such that x(3,i) < 0 corresponds
%     to a point that is farther from the objective than the focal plane.
%   e is a 3-by-n_rays matrix, where each column is a unit vector giving
%     the propagation direction
%   fobj is a 2-vector, [ffront fback], of focal lengths for the objective
%     lens. (See GENERATE_RAYS for help in specifying this quantity.)
%   ftube is a 2-vector, [ffront fback], of focal lengths for the tube
%     lens.
%   separation is the distance between the back principal planes of the
%     objective and tube lenses.
% and
%   xo is the set of points in the image plane for each ray (a
%     3-by-n_rays matrix). xo(3,i) = 0 defines the image plane.
%   eo is the matrix containing the output propagation directions.
%   s is a vector of optical path lengths traversed by each ray to get to
%     the output point xo.
%
% See also: TRACE_RAYS_SINEOBJ, GENERATE_RAYS.

% Copyright 2008 by Timothy E. Holy
  % First trace to focal plane
  dfp = -x(3,:) ./ e(3,:);
  x = x + repmat(dfp,3,1) .* e;
  % Trace to the conjugate plane
  M = -fb(2)/fa(2);
  xo = M*x;
  etilde = e(1:2,:);
  e1tilde = -x(1:2,:)/fa(2);
  eotilde = (-fa(1)/fb(1))*etilde + (1/fb(1)) * (x(1:2,:) - xo(1:2,:) + ...
    repmat(separation./sqrt(1-sum(e1tilde.^2,1)),2,1).*e1tilde); 
  eo = [eotilde; sqrt(1-sum(eotilde.^2,1))];
  if (nargout > 2)
    e1t3 = sqrt(1-sum(e1tilde.^2,1));
    s = dfp + separation./e1t3 - fa(1) + fb(2);
end
