function [n,mu,C] = cn_moments(x,covarianceModel,w)
% cn_moments: compute moments as the # of included points grows
%
% Syntax:
%   [n,mu,C] = cn_moments(x,covarianceModel)
%   [n,mu,C] = cn_moments(x,covarianceModel,w)
% where
%   x is a d-by-N matrix containing N data points (as columns) in d
%     dimensions, in the order you want them to be considered
%   covarianceModel is 'isotropic', 'diagonal', or 'full' (see
%     dist_mahalanobis)
%   w (optional) is a weight/multiplicity for each point
% and
%   n is a 1-by-N vector, giving to cumulative number of points (or the sum
%     of the weights, if supplied)
%   mu is d-by-N matrix storing the (weighted) mean of each neighborhood
%   C contains the sample covariance matrices, normalized by n-1 (or d*n-1
%     in the case of covarianceModel = 'isotropic'). The shape and size
%     depends upon the covarianceModel.
%
% See also: dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy

  [d,N] = size(x);
  if (nargin < 3)
    % Version with no weights
    n = 1:N;
    mu = bsxfun(@rdivide,cumsum(x,2),n);
    if (nargout > 2)
      switch covarianceModel
        case 'isotropic'
          x2 = sum(x.^2,1);
          mu2 = sum(mu.^2,1);
          C = (cumsum(x2) - n.*mu2)./(n*d-1);
        case 'diagonal'
          C = bsxfun(@rdivide,cumsum(x.^2,2) - bsxfun(@times,n,mu.^2),n-1);
        case 'full'
          xx = bsxfun(@times,reshape(x,[d 1 N]),reshape(x,[1 d N]));
          mumu = bsxfun(@times,reshape(mu,[d 1 N]),reshape(mu,[1 d N]));
          nrs = reshape(n,[1 1 N]);
          C = bsxfun(@rdivide,cumsum(xx,3) - bsxfun(@times,nrs,mumu),nrs-1);
        otherwise
          error('Covariance model not recognized');
      end
    end
  else
    error('Weights not yet implemented');
    % Version with weights
    n = sum(w);
    xw = bsxfun(@times,x,w);
    mu = sum(xw,2)/n;
    if (nargout > 2)
      dx = bsxfun(@minus,x,mu);
      dxw = bsxfun(@times,dx,w);
      switch covarianceModel
        case 'isotropic'
          C = sum(dx(:).*dxw(:))/(n*d-1);
        case 'diagonal'
          C = sum(dx.*dxw,2)/(n-1);
        case 'full'
          C = (dx*dxw')/(n-1);
        otherwise
          error('Covariance model not recognized');
      end
    end
  end
