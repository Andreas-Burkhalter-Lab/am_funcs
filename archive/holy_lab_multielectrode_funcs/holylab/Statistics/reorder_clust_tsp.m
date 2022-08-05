function [clusterOrder,minD] = reorder_clust_tsp(varargin)
% REORDER_CLUST_TSP: sort clusters by minimum "travel" distance
%
% This algorithm orders clusters according to the shortest "path" needed to
% sequentially "visit" each cluster. This is therefore a traveling salesman
% problem (tsp). This calls a genetic algorithm to find a good solution,
% but be aware that the solution may vary somewhat from call to call. You
% can call this repeatedly and monitor the total distance to choose the
% best solution.
%
% Syntax:
%   [clusterOrder,minD] = reorder_clust_tsp(D,options)
%   [clusterOrder,minD] = reorder_clust_tsp(X,clust,options)
% where
%  D is a matrix of distances between clusters
% OR
%  X is a d-by-N matrix of coordinates of raw data points
%  clust is a 1-by-N vector giving the cluster identity of these points
% and options may have the fields
%   pop_size (default 1000): see tspo_ga
%   num_iter (default 100): see tspo_ga
%   penaltyfunc: if present, should have the syntax
%       val = penaltyfunc(indx)
%     where indx is an ordering of 1:n_clusters, and val is a scalar. You
%     can use this to favor some orderings over others. For example,
%       penaltyfunc = @(indx) 20*(indx(1) ~= 5)
%     would penalize any ordering that doesn't start with cluster #5.
%     Alternatively supply this as a structure with two fields:
%        type: a string, currently supported is 'maxmean'
%        amplitude: a scalar
%     The 'maxmean' type is defined as
%        val = @(indx) amplitude*sum(abs(maxmeanX(indx) - sort(maxmeanX(indx),'descend')));
%     where maxmeanX is the maximum coordinate of the mean for each
%     cluster.  This penalty function favors putting clusters with the
%     largest mean coordinate values first. This is available only with the
%     (X,clust) syntax.
% and
%   clusterOrder is an ordering of 1:n_clusters that represents the
%     shortest path found; see reorder_clust to implement the change in
%     cluster #.
%   minD is the length of this shortest path.
%
% See also: reorder_clust, tspo_ga.

% Copyright 2010 by Timothy E. Holy

  options = struct;
  if (length(varargin) < 2)
    D = varargin{1};
  else
    if isnumeric(varargin{2})
      X = varargin{1};
      clust = varargin{2};
      [D,meanX] = rct_calcdist(X,clust);
    else
      D = varargin{1};
    end
  end
  if isstruct(varargin{end})
    options = varargin{end};
  end
  
  options = default(options,'pop_size',1000,'num_iter',100,'penaltyfunc',@(z) 0);
  
  if isstruct(options.penaltyfunc)
    switch options.penaltyfunc.type
      case 'maxmean'
        maxmeanX = max(meanX,[],1);
        options.penalty = @(indx) options.penaltyfunc.amplitude*sum(abs(maxmeanX(indx) - sort(maxmeanX(indx),'descend')));
      otherwise
        error(['Type ' options.penaltyfunc.type ' not recognized'])
    end
  end
  
  [clusterOrder,minD] = tspo_ga(zeros(size(D,1),0),D,options.pop_size,options.num_iter,false,false,options.penaltyfunc);

function [D,meanX] = rct_calcdist(X,clust)
  clabel = agglabel(clust);
  n_clust = length(clabel);
  meanX = zeros(size(X,1),n_clust);
  for i = 1:n_clust
    meanX(:,i) = mean(X(:,clabel{i}),2);
  end
  D = sqrt(sqrdist(meanX,meanX));