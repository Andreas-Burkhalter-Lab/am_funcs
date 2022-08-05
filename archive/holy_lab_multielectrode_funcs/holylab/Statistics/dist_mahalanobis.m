function z2 = dist_mahalanobis(dx,C,covarianceModel)
% dist_mahalanobis: a generalized distance based on a covariance
%
% Syntax:
%   z2 = dist_mahalanobis(dx,C,covarianceModel)
% Here, dx is a d-by-N matrix of displacements from the base point, for N
% points in d dimensions.
%
% This function permits several different parametrizations of constrained
% covariance matrices, as indicated by the optional flag covarianceModel:
%  'isotropic' is a covariance matrix proportional to the identity matrix,
%     C = sigma^2 * eye(d,d). In this case C is the scalar sigma^2.
%  'diagonal' is a covariance matrix with all off-diagonal elements zero,
%     but permits different values in each element of the diagonal. In this
%     case C is a column vector of length d.
%  'full' implies a full d-by-d covariance matrix.
%
% On output, z2 is the "squared z score", equal to
%   z2(i) = dx(:,i)'*(C\dx(:,i))
% although it is not computed that way.
%
% See also: moments_from_points.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 3)
    if isscalar(C)
      covarianceModel = 'isotropic';
    elseif isvector(C)
      covarianceModel = 'diagonal';
    else
      covarianceModel = 'full';
    end
  end
  
  switch covarianceModel
    case 'isotropic'
      z2 = sum(dx.^2,1)/C;
    case 'diagonal'
      z2 = sum(bsxfun(@rdivide,dx.^2,C),1);
    case 'full'
      R = chol(C);
      z2 = R'\dx;
      z2 = sum(z2.^2,1);
    otherwise
      error('Covariance model not recognized');
  end
