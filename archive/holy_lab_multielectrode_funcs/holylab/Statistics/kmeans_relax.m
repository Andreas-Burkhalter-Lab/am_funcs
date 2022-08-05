function [centroids,cindex,md] = kmeans_relax(X,centroids,cindexold)
% kmeans_relax: optimize the position of centroids, given a starting guess
%
% This function moves centroids using the k-means algorithm. You supply the
% initial positions, and this algorithm "flows" them to a (local) optimum.
%
% Syntax:
%   [centroidsOut,cindex,md] = kmeans_relax(X,centroidsIn)
%   [centroidsOut,cindex,md] = kmeans_relax(X,centroidsIn,cindexIn)
% where
%   X is the initial d-by-N data matrix, each d-dimensional point stored in
%     one column of X.
%   centroidsIn is the d-by-k matrix of intial centroid positions
% and
%   centroidsOut contains the final matrix of centroid positions
%   cindex is a 1-by-N vector, for each data point giving the index of the
%     nearest centroid
%   md is a 1-by-N vector, for each data point giving the distance to the
%     nearest centroid.
%
% The optional input cindexIn can be a previous centroid assignment vector;
% if you are doing some form of iterative or incremental k-means, supplying
% this can  potentially save a small amount of time, because it allows the
% algorithm to skip updating any centroids whose position should not
% change. This input affects only the first iteration of k-means, so in
% many circumstances supplying it will have negligible effect.
%
% See also: em_relax.

% Copyright 2010 by Timothy E. Holy

  %% Initialization
  if (nargin > 1)
    n_centroids = size(centroids,2);
  else
    n_centroids = 0;
  end
  if (n_centroids < 2)
    % 0 or 1 centroids: just return the mean of the data set
    centroids = mean(X,2);
    cindex = ones(1,size(X,2));
    if (nargout > 2)
      md = sum(bsxfun(@minus,X,centroids).^2,1);
    end
    return
  end
  % Find closest centroid
  [md,cindex] = mindist(X,centroids); % compute closest centroid & distance
  % Set the "moved" flag, which indicates which centroids have had their
  % membership changed by at least one data point
  if (nargin < 3 || isempty(cindexold))
    moved = true(1,n_centroids);     % all centroids need updating on the first iteration
  else
    moved = kmr_update_moved(cindex,cindexold,n_centroids);
  end
  
  %% Run the k-means algorithm
  while any(moved)
    clabel = agglabel(cindex);  % separate points into groups assigned to a given centroid
    for i = 1:length(clabel)
      % Note we don't want to assume that all centroids still have points
      % assigned to them
      if (moved(i) && ~isempty(clabel{i}))
        centroids(:,i) = mean(X(:,clabel{i}),2);  % Recompute centroid position
      end
    end
    cindexold = cindex;
    [md,cindex] = mindist(X,centroids);  % Reassign points
    moved = kmr_update_moved(cindex,cindexold,n_centroids);
  end
  
  %% Remove any empty centroids
  l = cellfun(@length,clabel);
  keepFlag = l > 0;
  centroids = centroids(:,keepFlag);
  cskf = cumsum(keepFlag);
  cindex = cskf(cindex);
  
function moved = kmr_update_moved(cindex,cindexold,n_centroids)
  moved = false(1,n_centroids);
  ischanged = cindexold ~= cindex;
  moved(cindexold(ischanged)) = true;
  moved(cindex(ischanged)) = true;
    