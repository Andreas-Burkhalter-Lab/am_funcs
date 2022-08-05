function [A,mu,sigma,fx] = zernproj2gaussian(p,normalizedFlag)
% zernproj2gaussian: convert projections onto Zernike basis into gaussian parameters
%
% Suppose we have an aberration characterized by its projection onto
% Zernike functions of first and second order. We can ask what
% Gaussian-shaped hump best accounts for these projections. Because an
% isotropic 2-d gaussian has 4 parameters (amplitude, two for the centroid,
% and one for the standard deviation) and there are 5 projections, the
% system of equations is potentially overdetermined (which can be useful for
% handling noise).
%
% Syntax:
%   [A,mu,sigma] = zernproj2gaussian(p,normalizedFlag)
% where
%   p is a 5-vector containing the projections onto Z(1,-1), Z(1,1),
%     Z(2,-2), Z(2,0), and Z(2,2), repsectively;
%   normalizedFlag (default true) specifies whether the projections were
%     calculated using normalized Zernikes (see ZERNIKE_VALUES)
% and
%   A is the amplitude of the normalized gaussian;
%   mu is a 2-vector containing the position of the centroid;
%   sigma is the estimate of the width of the gaussian.
%
% To be specific, the model is fitting to the functional form
%      g(x,y) = A / (2*pi*sigma^2) * exp(-((x-mux)^2 + (y-muy)^2)/(2*sigma^2))
%
% See also: ZERNIKE_VALUES.

% Copyright 2009 by Timothy E. Holy

  %% Input parsing
  if (nargin < 3)
    normalizedFlag = true;
  end
  if (length(p) ~= 5)
    error('p must be a vector of length 5');
  end
  if normalizedFlag
    p = p(:) ./ [2; 2; sqrt(6); sqrt(3); sqrt(6)];
  end
  
  %% Initial guess
  % Assuming non-normalized Zernikes, the projections of this Gaussian are:
  %   A*mux
  %   A*muy
  %   A*(mux^2 - muy^2)
  %   A*(4*sigma^2 + 2*mux^2 + 2*muy^2 - 1)
  %   A*2*mux*muy
  % (Be aware that 3 and 5 can be swapped depending on the convention (x,y)
  % or (y,x) for matrices<->2d functions)
  % We want to solve these equations for A, mux, muy, and sigma in a
  % least-squares sense.
  % 4 of the equations are independent of sigma. Solve these first, then
  % get sigma.
  
  if (abs(p(5)) > abs(p(3)))
    Aguess = 2*p(1)*p(2)/p(5);
  else
    Aguess = (p(1)^2 - p(2)^2)/p(3);
  end
  muguess = p(1:2)/Aguess;
  x0 = [Aguess; muguess(:)];
  x = fsolve(@(x) projmatch(x,p),x0,optimset('Jacobian','on','NonlEqnAlgorithm','gn'));
  A = x(1);
  mu = x(2:3);
  fx = projmatch(x,p);
  sigma2 = (p(4)/A - 2*sum(mu.^2) + 1)/4;
  if (sigma2 < 0)
    warning('phasediversity:zernproj2gaussian','negative value for sigma found, truncating to 0');
    sigma2 = 0;
  end
  sigma = sqrt(sigma2);
end

function [f,J] = projmatch(x,p)
  A = x(1);
  mux = x(2);
  muy = x(3);
  f(1) = A*mux - p(1);
  f(2) = A*muy - p(2);
  f(3) = A*(mux^2 - muy^2) - p(3);
  f(4) = A*2*mux*muy - p(5);
  if (nargout > 1)
    J = zeros(4,3);
    J(1,:) = [mux, A, 0];
    J(2,:) = [muy, 0, A];
    J(3,:) = [(mux^2 - muy^2), 2*A*mux, 2*A*muy];
    J(4,:) = [2*mux*muy, A*2*muy, A*2*mux];
  end
end
