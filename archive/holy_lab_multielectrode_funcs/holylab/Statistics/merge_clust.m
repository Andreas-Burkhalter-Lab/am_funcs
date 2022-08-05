function [clust,clabel,nlabel] = merge_clust(clabel,merge_groups)
% MERGE_CLUST: manually merge clusters by their cluster #
% Syntax:
%   clust = merge_clust(clustold,merge_groups)
%   clust = merge_clust(clabel,merge_groups)
%   [clust,clabel,nlabel] = merge_clust(...)
% where
%   clustold contains the input cluster #s
% OR
%   clabel is a cell array, each entry specifying the indices of the points
%     in a given cluster;
% and
%   merge_groups is a cell array, each entry specifying a list of clusters
%     you want to combine.
% and
%   clust is a new cluster label, one for each point
%   clabel is the new output clabel after merging
%   nlabel is the # of points in each new cluster
%
% For example, if you have 23 clusters, and you want to merge clusters 3
% and 17 together, and clusters 5, 6, and 8 together, then you would do
% this:
%   [clustnew,clabelnew,nlabelnew] = merge_clust(clabel,{[3 17],[5 6 8]})
% Note that the clusters are numbered in a way that is consistent with the
% original ordering, and all merge groups are combined in with the first
% one on the list.
%
% See also: AGGLABEL.

% Copyright 2010 by Timothy E. Holy

  n_merge = length(merge_groups);
  if isnumeric(clabel)
    % clust input
    clabel = agglabel(clabel);
  end
  
  % For each merge group, combine all into the first one in the list
  for i = 1:n_merge
    thisMerge = merge_groups{i};
    clabel_merge = cat(2,clabel{thisMerge});
    clabel{thisMerge(1)} = clabel_merge;
    for j = thisMerge(2:end)
      clabel{j} = [];
    end
  end
  % Eliminate empty groups
  l = cellfun(@length,clabel);
  clabel = clabel(l > 0);
  nlabel = l(l > 0);
  % Reassign cluster numbers
  clust = zeros(1,sum(nlabel));
  for i = 1:length(clabel)
    clust(clabel{i}) = i;
  end
      