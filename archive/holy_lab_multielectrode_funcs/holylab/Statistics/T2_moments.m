function [mu1,mu2,mu3] = T2_moments(d,n,covtype)
% T2_moments: compute the moments of T^2 with constrained covariance matrices
% Given data points drawn from a d-dimensional Gaussian centered at zero,
% this computes the first three central moments of
%   T^2 = n xbar'*inv(W)*xbar,
% where
%   n is the number of data points drawn from the Gaussian, and W is the
%   "covariance matrix." If we let X be the d-by-n data matrix, and define
%    Xc as the centered d-by-n data matrix (conceptually Xc = X -
%    mean(X,2)), then W may be the full sample covariance,
%       W = Xc*Xc'/(n-1)
%   or a diagonal approximation,
%       W = diag(sum(Xc.^2,2))/(n-1)
%   or the isotropic variant proportional to the identity,
%       W = eye(d,d)*sum(Xc(:).^2)/(d*n-1).
%
% Syntax:
%   [mu1,mu2,mu3] = T2_moments(d,n,covtype)
% where
%    d is the dimensionality
%    n is the number of data points
%    covtype = 0 is the isotropic covariance, 1 is the diagonal covariance,
%      and 2 is the full covariance
% and
%    mu1 is the expected value of T^2
%    mu2 is the expected variance of T^2
%    mu3 is the expected third centered moment, E((T^2-mu1)^3).
%
% Either d or n may be a vector, but not both.
%
% Note that the moments are not defined for certain parameter combinations,
% and if they are defined, they are not guaranteed to be finite. These
% problems disappear for larger n.
%
% See also: fdistribution_approx_params.

% Copyright 2011 by Timothy E. Holy

  switch covtype
    case 0
      if (nargout > 2 && d*(n-1) < 6)
        error('Third moment not defined');
      elseif (nargout > 2 && d*(n-1) < 4)
        error('Second moment not defined');
      elseif (d*(n-1) < 2)
        error('First moment not defined');
      end
      mu1 = d.*(d*n-1)./(d*(n-1)-2);
      mu2 = 2*d.*(d*n-1).^2.*(d*n-2)./(d*(n-1)-4)./(d*(n-1)-2).^2;
      mu3 = 8*d.*(d*n-1).^3.*(d*n-2).*(d*(n+1)-2)./(d*(n-1)-6)./(d*(n-1)-4)./(d*(n-1)-2).^3;
    case 1
      if (nargout > 2 && n < 7)
        error('Third moment not defined');
      elseif (nargout > 2 && n < 5)
        error('Second moment not defined');
      elseif (n < 3)
        error('First moment not defined');
      end
      mu1 = d*(n-1)./(n-3);
      mu2 = 2*d*(n-1).^2.*(n-2)./(n-5)./(n-3).^2;
      mu3 = 8*d*(n-2).*(n-1).^4./(n-7)./(n-5)./(n-3).^3;
    case 2
      if (nargout > 2 && n < d+6)
        error('Third moment not defined');
      elseif (nargout > 2 && n < d+4)
        error('Second moment not defined');
      elseif (n < d+2)
        error('First moment not defined');
      end
      mu1 = d*(n-1)./(n-d-2);
      mu2 = 2*d*(n-1).^2.*(n-2)./(n-d-4)./(n-d-2).^2;
      mu3 = 8*d*(n-2).*(n-1).^3.*(n+d-2)./(n-d-6)./(n-d-4)./(n-d-2).^3;
    otherwise
      error('covtype unkown');
  end
  