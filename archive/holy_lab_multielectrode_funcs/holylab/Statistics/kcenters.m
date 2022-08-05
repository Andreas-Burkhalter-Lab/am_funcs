function [indx,c] = kcenters(X,k)
% KCENTERS: the k-centers/k-medoids/PAM algorithm
% Syntax:
%   [indx,c] = kcenters(X,k)
% where
%   X is the data matrix (d-by-N, d = # of variables, N = # observations)
%   k is the desired number of clusters
% and
%   indx is the cluster index
%   c is the d-by-k matrix containing the medoids
  
% This is just a wrapper around Kmedoid
  
  data.X = X';
  param.c = k;
  param.vis = 0;

  r = Kmedoid(data,param);
  
  [mx,indx] = max(r.data.f,[],2);
  indx = indx';
  c = r.cluster.v';
  