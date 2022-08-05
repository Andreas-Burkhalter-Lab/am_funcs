function [err,x0] = ray_convergence_error(xtot,stot)
% RAY_CONVERGENCE_ERROR: compute aberrations given a set of rays
%
% For a set of rays, where ray i passes through position x(i) (a vector)
% and propagates in the direction s(i) (another vector), the convergence
% error is defined as the minimum with respect to x0 of the function
%
%     err = mean_i [(x(i)-x0)^2  - ( (x(i)-x0) . s(i) )^2]
%
% where the '.' means dot product and x0 is a fixed point. In other
% words, this measures the positional error assuming the rays pass
% through a single point x0.  The point x0 is chosen to minimize this
% error.
%
% Syntax:
%   [err,x0] = ray_convergence_error(r)
% where r is a structure array of rays (see RAY)
% OR
%   [err,x0] = ray_convergence_error(x,s)
% where x is a d-by-N matrix, each column giving one point on each ray,
% and s is a d-by-N matrix, each column giving the propagation direction
% of each ray. (d = 2 or 3, depending on the dimensionality of the optics
% problem).
%
% See also: RAY.

% Copyright 2006-2008 by Timothy E. Holy
  
  if (nargin == 1)
    r = xtot([xtot.valid]);  % restrict to valid rays
    xtot = [r.x0];
    stot = [r.e];  % the ray structure uses "e" for "s" (confusing, I know)
  end
  [n_dims,n_rays] = size(xtot);
  if (n_rays < 2)
    err = 0;
    x0 = nan(size(xtot));
    return
  end
  % Set up the SVD solution
  b = zeros(n_dims*n_rays,1);
  A = zeros(n_dims*n_rays,n_dims);
  I = eye(n_dims,n_dims);
  for k = 1:n_rays
    indx = (k-1)*n_dims + (1:n_dims);
    x = xtot(:,k);
    e = stot(:,k);
    xdote = x' * e;
    b(indx) = x - xdote * e;
    A(indx,:) = I - e*e';
  end
  % Calculate the point of convergence
  x0 = A\b;
  % Test to see if any rays are imaginary
  if any(b ~= real(b))
    err = inf;   % Give infinite error
  else
    dx = b - A*x0;
    err = sum(dx.^2)/n_rays;
  end
