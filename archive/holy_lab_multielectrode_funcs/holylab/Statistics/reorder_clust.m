function [clust,clabel,nlabel] = reorder_clust(clust,clusterOrder,clabel,nlabel)
% REORDER_CLUST: re-label clusters, arranging them in a different order
% Syntax:
%   clust = reorder_clust(clust,clusterOrder)
%   [clust,clabel,nlabel] = reorder_clust(clust,clusterOrder,clabel,nlabel)
% where
%   clust is the cluster number assigned to each point
%   clusterOrder is the new order in which you want clusters to appear
%   clabel, nlabel are the outputs of agglabel
% and the outputs are like the inputs, but with new cluster numbers
% assigned.
%
% clusterOrder might, of course, be assigned by an automated algorithm,
% such as reorder_clust_tsp.
%
% See also: AGGLABEL, REORDER_CLUST_TSP.

% Copyright 2010 by Timothy E. Holy

if (length(unique(clusterOrder)) ~= max(clust(:)))
  error('clusterOrder must have one entry for each cluster in clust');
end
[~,clusterOrderInv] = sort(clusterOrder);
clust = clusterOrderInv(clust);
if (nargin > 2)
  clabel = clabel(clusterOrder);
end
if (nargin > 3)
  nlabel = nlabel(clusterOrder);
end
