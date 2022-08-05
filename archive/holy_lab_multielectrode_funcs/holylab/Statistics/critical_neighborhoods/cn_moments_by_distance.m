function [sortOrder,n,mu,C] = cn_moments_by_distance(x,varargin)
% Like cn_moments except ordering is determined by increasing distance. But
% note the mean is computed relative to the basepoint x0.
% Syntax:
%   [sortOrder,n,mu,C] = cn_moments_by_distance(x,params,x0)
%   [sortOrder,n,mu,C] = cn_moments_by_distance(x,params,x0,C0)
%   [sortOrder,n,mu,C] = cn_moments_by_distance(x,w,params,...)
% w is a weight vector. C0 needs to be supplied only for
%    params.covarianceModel ~= 'isotropic'.

  indx = 1;
  useWeight = false;
  if ~isstruct(varargin{1})
    w = varargin{1};
    indx = 2;
    useWeight = true;
  end
  params = varargin{indx};
  x0 = varargin{indx+1};
  if length(varargin) > indx+1
    C0 = varargin{indx+2};
  end
  
  %% Order the points
  dx = bsxfun(@minus,x,x0);
  if strcmp(params.covarianceModel,'isotropic')
    [~,sortOrder] = sort(sum(dx.^2,1));
  else
    dist = dist_mahalanobis(dx,C0,params.covarianceModel);
    [~,sortOrder] = sort(dist);
  end

  %% Calculate the moments
  if useWeight
    [n,mu,C] = cn_moments(dx(:,sortOrder),params.covarianceModel,w);
  else
    [n,mu,C] = cn_moments(dx(:,sortOrder),params.covarianceModel);
  end
end
