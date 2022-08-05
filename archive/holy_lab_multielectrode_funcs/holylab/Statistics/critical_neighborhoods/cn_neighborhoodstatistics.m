function T2thresh = cn_neighborhoodstatistics(N,d,pvalue,covarianceModel,nMin)
% cn_neighborhoodstatistics: calculate threshold for significance of T^2
%
% Given a set of points, T^2 is a statistic for the displacement of their
% mean.  This function calculates the threshold for statistical
% significance for neighborhoods of all sizes (i.e., number of points) up
% to a specified maximum.
%
% Syntax:
%   T2thresh = cn_neighborhoodstatistics(N,d,pvalue,covarianceModel)
%   T2thresh = cn_neighborhoodstatistics(N,d,pvalue,covarianceModel,nMin)
% where
%   N is the maximum number of points
%   d is the dimensionality
%   pvalue is the threshold for significance, e.g., 0.01
%   covarianceModel is a flag indicating the parametrization of the
%     covariance, see dist_mahalanobis.
%   nMin (optional) is the minimum number of points in a neighborhood;
%     the first nMin values of T2thresh are set to infinity
% and
%   T2thresh is a vector containing the threshold for significance for
%     neighborhoods with 1:N points.
%
% See also: dist_mahalanobis, cn_preordered.
  
% Copyright 2011 by Timothy E. Holy
  
  n = 1:N;
  isbad = false(1,N);
  switch covarianceModel
    case 'isotropic'
      n1 = d;
      n2 = d*(n-1);
      s = (n-1)./(d*n-1);
      isbad(1) = true;
    case 'diagonal'
      % Use the 2-moment approximation
      n1 = d;
      n2 = ((d+2)*n-5*d+2)/3;
      s = (n-3).*n2./(d*(n-1).*(n2-2));
      isbad(1:4) = true;
    case 'full'
      n1 = d;
      n2 = n-d;
      s = (n-d)./(d*(n-1));
      isbad(1:d) = true;
    otherwise
      error('Covariance model not recognized')
  end
  T2thresh = finv(1-pvalue,n1,n2)./s;
  T2thresh(isbad) = inf;
  if (nargin > 4)
    T2thresh(1:nMin) = inf;
  end
  
