function [n,mu,C] = moments_from_points(x,covarianceModel,w)
% moments_from_points: compute zeroth through 2nd moments from a set of points
%
% Syntax:
%   [n,mu,C] = moments_from_points(x,covarianceModel)
%   [n,mu,C] = moments_from_points(x,covarianceModel,w)
% where
%   x is a d-by-N matrix containing N data points (as columns) in d
%     dimensions
%   covarianceModel is 'isotropic', 'diagonal', or 'full' (see
%     dist_mahalanobis)
%   w (optional) is a weight/multiplicity to be supplied for each point
% and
%   n is the number of points, or the sum of the weights if supplied
%   mu is the mean of the (weighted) points
%   C is the sample covariance matrix, normalized by n-1 (or d*n-1 in the
%     case of covarianceModel = 'isotropic')
%
% See also: dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy

  [d,n] = size(x);
  if (nargin < 3)
    % Version with no weights
    mu = mean(x,2);
    if (nargout > 2)
      dx = bsxfun(@minus,x,mu);
      switch covarianceModel
        case 'isotropic'
          C = sum(dx(:).^2)/(n*d-1);
        case 'diagonal'
          C = sum(dx.^2,2)/(n-1);
        case 'full'
          C = (dx*dx')/(n-1);
        otherwise
          error('Covariance model not recognized');
      end
    end
  else
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
