function [b,gXb,gRb] = biharmonic_pointsource(varargin)
% BIHARMONIC_POINTSOURCE: shape of membrane under point force
%
% This function computes the shape of a membrane that obeys the
% biharmonic equation on a circular boundary with Dirichlet boundary
% conditions, i.e.,
%
%     laplacian^2 b = 0,  with b = 0 and db/dn = 0 on the circle
%
% One can find the relevant Green's function quoted in Y. A. Melnikov,
% "Influence Functions of a Point Force for Kirchhoff Plates with Rigid
% Inclusions," Journal of Mechanics, 20, pp. 249-256 (2004).  This paper
% in turn gives credit to a reference I haven't pursued, N. N. Lebedev,
% I. M. Skal'skaya, and Ya. S. Uflyand, Problems of Mathematical Physics,
% Pergamon Press, New York (1966).
%
% Syntax:
%   b = biharmonic_pointsource(n,x0)
% n is a scalar (# of pixels along each axis), x0 is a 1-by-2 vector
% containing the position of the point source, and we assume that the
% radius R of the circle is 1. (Note sum(x0.^2) must be less than R).
%
%   b = biharmonic_pointsource(X,x0,R)
% X is a N-by-2 matrix of points at which to evaluate the surface
%
%   [b,gXb,gRb] = biharmonic_pointsource(X,x0,R)
% also computes the gradient of b with respect to X and R.
  
% Copyright 2009 by Timothy E. Holy
  
  x0 = varargin{2};
  if (isscalar(varargin{1}))
    % "Image" mode
    x0 = x0([2 1]);  % because x & y are flipped in images
    x0(1) = -x0(1);  % because ydir is 'reverse'
    n = varargin{1};
    b = zeros(n,n);
    x = (1:n) - (n+1)/2;
    x = x/max(abs(x));
    [X{1},X{2}] = ndgrid(x,x);
    R2 = X{1}.^2 + X{2}.^2;
    pupil = R2 < 1;
    Xp = cell(1,2);
    dXp = cell(1,2);
    for dimIndex = 1:2
      Xp{dimIndex} = X{dimIndex}(pupil);
      dXp{dimIndex} = Xp{dimIndex} - x0(dimIndex);
    end
    ct{1} = 1 - (Xp{1}*x0(1) + Xp{2}*x0(2));
    ct{2} = - (Xp{2}*x0(1) - Xp{1}*x0(2));
    dx2 = dXp{1}.^2 + dXp{2}.^2;
    ct2 = ct{1}.^2 + ct{2}.^2;
    lt = log(dx2 ./ ct2);
    lt(dx2 == 0) = 0;
    b(pupil) = dx2 .* lt + (1-R2(pupil)) * (1 - sum(x0.^2));
  else
    % "Coordinate" mode
    % In the literature the formula is presented in terms of complex
    % numbers, but here we implement everything in terms of real numbers
    R = varargin{3};
    X = varargin{1}/R;
    x0 = x0/R;
    X2 = sum(X.^2,2);
    x02 = sum(x0.^2);
    inpupil = X2 < 1;
    b = zeros(size(X,1),1);
    X = X(inpupil,:);
    X2 = X2(inpupil);
    %N = size(X,1);
    dX = [X(:,1)-x0(1),X(:,2)-x0(2)];
    Xdotx0 = X(:,1)*x0(1) + X(:,2)*x0(2);
    Xcrossx0 = X(:,2)*x0(1) - X(:,1)*x0(2);  % cross product (z-component)
    dX2 = sum(dX.^2,2);
    ct2 = (1-Xdotx0).^2 + Xcrossx0.^2;
    lt = log(dX2 ./ ct2);
    lt(dX2 == 0) = 0;
    b(inpupil) = dX2.*lt + (1-X2)*(1-x02);
    %b(inpupil) = (1-X2)*(1-x02);
    if (nargout > 1)
      % Gradient with respect to X
      gXb = zeros(size(b,1),2);
      coef = 2*dX2./ct2;
      coefa = coef.*(1-Xdotx0);
      coefb = coef.*Xcrossx0;
      gXb(inpupil,:) = 2*dX.*repmat(1+lt,1,2) - 2*(1-x02)*X + ...
        [coefa*x0(1),coefa*x0(2)] + [coefb*x0(2),-coefb*x0(1)];
      gXb = gXb/R;
      % Gradient with respect to R
      gRb = zeros(size(b));
      gRb(inpupil) = (2/R)*(dX2.*(1-lt-2*(1-Xdotx0)./ct2)-2*(1-X2)*(1-x02)+2-X2-x02);
    end
  end
