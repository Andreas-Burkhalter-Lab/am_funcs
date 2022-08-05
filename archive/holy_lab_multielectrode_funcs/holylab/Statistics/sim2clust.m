function clust = sim2clust(Y,thresh)
% SIM2CLUST: convert similarity matrix into a set of clusters
% Syntax:
%   clust = sim2clust(Y,thresh)
% where
%   Y is the type of ouput produced by MAPS2SIM
%   thresh is a threshold for joining classes
% and
%   clust is a vector indicating the cluster number assigned to each
%     row(/column) of Y.
%
% See also: MAPS2SIM.
  
% Copyright 2008 by Timothy E. Holy
  
  Y = (Y >= thresh);
  todo = true(1,size(Y,1));
  cindx = 1;
  clustc = {};
  while ~isempty(cindx)
    clustc_tmp = find(Y(cindx,:));
    clustc{end+1} = clustc_tmp;
    todo(clustc_tmp) = false;
    cindx = find(todo,1,'first');
  end
  clust = zeros(size(todo));
  for i = 1:length(clustc)
    clust(clustc{i}) = i;
  end
    