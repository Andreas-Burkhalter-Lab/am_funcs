function [v]=proj_with_svd(x)
   n_pts = size(x,2);
   xbar = mean(x,2);
   dx = x - repmat(xbar,[1 n_pts]);
   [U,S,V] = svd(dx,'econ');
   proj = U(:,1:n_pts-1);
   v=proj;
