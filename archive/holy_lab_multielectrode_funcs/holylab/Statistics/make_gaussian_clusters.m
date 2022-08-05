function [x,clust] = make_gaussian_clusters(n,xc,sigma)
% MAKE_GAUSSIAN_CLUSTERS: create simulated data containing clusters
% Syntax:
%   [x,clust] = make_gaussian_clusters(n,xc,sigma)
% where
%   n is a 1-by-n_clusters vector specifying the number of data points in
%     each cluster
%   xc is a d-by-n_clusters vector specifying the cluster center
%   sigma can be:
%     a scalar, in which case it is the (isotropic) standard deviation for
%       all gaussians
%     a vector of size 1-by-n_clusters, in which case each gaussian is
%       isotropic but may have different standard deviation
%     a tensor of size d-by-d-by-n_clusters, in which case sigma(:,:,i) is
%       the "square root of the covariance matrix" for the ith gaussian,
%       i.e., V*sqrt(D) where [V,D] = eig(C) (C is the covariance matrix).
% and
%   x is d-by-N, where N is the number of data points (=sum(n))
%   clust is 1-by-N, containing the cluster number associated with each
%     point.

% Copyright 2010 by Timothy E. Holy

n_clusters = length(n);
ncum = [0 cumsum(n)];
d = size(xc,1);
x = zeros(d,ncum(end));
clust = zeros(1,ncum(end));
for i = 1:n_clusters
  xr = randn(d,n(i));
  if isscalar(sigma)
    xr = sigma*xr;
  elseif isvector(sigma)
    xr = sigma(i)*xr;
  else
    xr = sigma(:,:,i)*xr;
  end
  xnew = bsxfun(@plus,xc(:,i),xr);
  rng = ncum(i)+1:ncum(i+1);
  x(:,rng) = xnew;
  clust(rng) = i;
end
