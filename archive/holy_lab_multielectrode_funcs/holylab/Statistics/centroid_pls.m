function [mu,iter] = centroid_pls(X,nu,order)
% CENTROID_PLS: find centroid so that power-law scaled data have zero mean
%
% Given a set of data vectors x_i and a scaling parameter nu, find the
% centroid mu such that the rescaled data vectors
%     y_i = (x_i - mu)/|x_i - mu|^(1-nu)
% have zero mean.
%
% Syntax:
%   mu = centroid_pls(X)
% where the x_i are the rows of X, assumes nu = 1 (no scaling, so
% mu = mean(X)); 
%   mu = centroid_pls(X,nu)
% allows you to explicitly specify nu;
%   mu = centroid_pls(X,nu,order)
% allows you to control the convergence of the algorithm, where order = 1
% results in a linearly-convergent algorithm with storage requirements of
% order Nd (N = # points, d = dimensionality), while order = 2 results in
% a quadratically-convergent algorithm with storage of order Nd + d^2.
% By default, order = 1, because the extra overhead of computing the
% Jacobian can overwhelm the convergence advantages of the quadratic
% algorithm, in addition to the possible storage issues.
%
% Finally,
%   [mu,iter] = centroid_pls(...)
% allows you to find out how many steps were required for convergence.
%
% See also: PCA.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    nu = 1;
  end
  if (nu == 1)
    mu = mean(X);
    iter = 0;
    return;
  end
  if (nargin < 3)
    order = 1;  % More steps to converge, but each step is faster
  end
  if (order ~= 1 && order ~= 2)
    error('order must be 1 or 2');
  end
  [N,d] = size(X);
  if (d > 100 && order == 2)
    warning('In large dimensions, may want order = 1');
  end
  % Initial value for mu
  mu = mean(X,1);
  xvar = mean(var(X,0,1));
  epsilon2 = xvar;  % We'll use epsilon to regularize denominators for
                   % points near mu
  iter = 1;
  while (epsilon2 > d*eps*xvar)
    dX = X - repmat(mu,N,1);
    dX2 = sum(dX.^2,2);
    nzIndex = find(dX2);
    % Compute RHS
    dXsc = dX;
    dXsc(nzIndex,:) = dXsc(nzIndex,:)./repmat(dX2(nzIndex).^((1-nu)/2),1,d);
    b = sum(dXsc,1);
    % Compute LHS
    denom = 1./(dX2 + epsilon2).^((1-nu)/2);
    if (order == 2)
      dXnorm = zeros(size(X));
      dXnorm(nzIndex,:) = dX(nzIndex,:)./repmat(sqrt(dX2(nzIndex)),1,d);
      A = sum(denom)*eye(d) - dXnorm'*(repmat((1-nu)*denom,1,d).*dXnorm);
    else
      A = sum(denom);
    end
    dmu = b/A;
    mu = mu + dmu;
    epsilon2 = sum(dmu.^2);
    iter = iter+1;
  end
