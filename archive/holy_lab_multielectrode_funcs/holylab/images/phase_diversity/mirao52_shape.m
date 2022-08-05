function [z,dzdv] = mirao52(X,v,s)
% MIRAO52: compute mirror shape from voltage
% Syntax:
%   z = mirao52(X,v,s)
%   [z,dzdv] = mirao52(X,v,s)
% where
%   X: each row of X is a different point in the pupil plane. You should
%     express these in _normalized coordinates_, e.g. rho .* cos(theta) and
%     rho .* sin(theta), so that the parameters of s will not depend upon
%     the vagaries of your situation (e.g., how many pixels are in the
%     region you're currently processing...).
%   v: a vector of voltages supplied to the 52 actuators
%   s has the following fields:
%     A: 3-by-2 matrix parameters of the affine transformation to apply to X
%        before computing z,
%            Xp = tformfwd(X,T) where T = maketform('affine',A)
%        (this includes all scaling, rotation, shearing, stretching)
%     sigma: width of region-of-influence of each actuator (gaussian)
%     v0 & v2z: given a voltage v(i) for the ith actuator, the gaussian is
%              z_i = (v(i)-v0(i))*v2z(i)*exp(-(xp-x0_i)^2/(2*sigma^2))
%       so v0 and v2z are vectors of length 52 OR are scalars
% and
%   z (a column vector) is the displacement at the corresponding position
%     of X.
%   dzdv is a npts-by-52 matrix containing the gradient of z with respect
%     to each actuator voltage.
%
% If you want to use this to compute the phase for all the points in the
% pupil: assume H0 is a pupil function (1 inside the pupil and 0
% outside), then do the following:
%   indx = find(H0(:));
%   [i,j] = ind2sub(size(H0),indx);
%   z = mirao52([i j],v,s);
%   phi = zeros(size(H0));
%   phi(indx) = (2*pi/lambda)*z;
  
% Copyright 2009 by Timothy E. Holy
  
  % This matrix specifies the geometry of the actuator grid
  act_loc = [0 0 11 19 27 35 0 0; ...
	     0 5 12 20 28 36 43 0; ...
	     1 6 13 21 29 37 44 49; ...
	     2 7 14 22 30 38 45 50; ...
	     3 8 15 23 31 39 46 51;...
	     4 9 16 24 32 40 47 52;...
	     0 10 17 25 33 41 48 0; ...
	     0 0 18 26 34 42 0 0];
  % Convert it to a matrix of x,y positions
  X0 = zeros(max(act_loc(:)),2);
  indx = find(act_loc);
  [i,j] = ind2sub(size(act_loc),indx);
  X0(act_loc(indx),:) = [i j];

  [n_points,d] = size(X);
  if (d ~= 2)
    error('Only supports 2-dimensional pupil functions')
  end
  calc_derivative = false;
  if (nargout > 1)
    dzdv = zeros(n_points,n_actuators);
    calc_derivative = true;
  end
  
  % Put the input pupil coordinates on the grid of the actuators
  T = maketform('affine',s.A);
  Xp = tformfwd(T,X);

  % Compute the coefficients
  z0 = (v(:) - s.v0(:)) .* s.v2z(:);
  
  % Calculate the displacement at each coordinate, 
  z = zeros(n_points,1);
  s2 = 2*s.sigma^2;
  for i = 1:n_actuators
    dX = Xp - repmat(X0(i,:),n_points,1);
    dX2 = sum(dX.^2,2);
    zia = exp(-dX2 / s2);
    z = z + z0(i)*zia;
    if calc_derivative
      dzdv(:,i) = s.v2z(i) * zia;
    end
  end
end