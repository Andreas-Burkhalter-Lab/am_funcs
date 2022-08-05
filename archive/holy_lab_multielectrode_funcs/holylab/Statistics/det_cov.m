function dC = det_cov(C,d,covarianceModel)
% det_cov: determinant of the covariance matrix, for different constrained parametrizations
%
% Syntax:
%   dC = det_cov(C,d,covarianceModel)
% where
%   C is the covariance matrix (see dist_mahalanobis for the format)
%   d is the dimensionality (needs to be supplied for the 'isotropic' case)
%   covarianceModel is 'isotropic', 'diagonal', or 'full'
% and
%   dC is the determinant of C.
%
% See also: dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy

  switch covarianceModel
    case 'isotropic'
      dC = C^d;
    case 'diagonal'
      dC = prod(C);
    case 'full'
      dC = det(C);
    otherwise
      error('Covariance model not recognized');
  end
end