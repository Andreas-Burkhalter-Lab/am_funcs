function [clust,results,options] = bn_cluster_slow(X,varargin)
% BN_CLUSTER_SLOW:  clustering using balanced neighborhood criterion
%
% Note: this was banged out quickly, and is intended as a temporary measure
% until faster code is ready. It is suitable only for relatively small
% problems (e.g., fewer than 10^4 points).
%
% Syntax:
%   clust = bn_cluster_slow(X)
%   clustY = bn_cluster_slow(X,Y)
%   [clust,results,options] = bn_cluster_slow(...,options)
% where
%   X is a d-by-N matrix of N data points in d dimensions
%   Y is a d-by-M matrix of M probe points in d dimension (if absent, Y
%     is set to X)
%   options is a structure which may have the following fields:
%     consolidate (default true if Y = X): if set to true, checks the
%       neighborhood of each "peak" to see if a plurality of neighbors
%       belong to a different group; if so, it reassigns the group.
%     n_iter_max (default 100): the maximum number of movement steps
%       allowed before deciding that the given point maps to itself
% and
%   clust is an integer index assigning the cluster number of each of the N
%     points.
%   results is a structure which has the following fields:
%     n_iter is the number of iterations used for each data point.
%     map is a 1-by-N vector that contains, for each point, the index of
%       the point at the "peak" of the uphill flow.
%     Ymap is a d-by-M vector that contains the locations of the probe
%       points at termination (note this is not flowed to the peak, only
%       until they are closer to another probe point's starting location
%       than its own)
%   options (on output) records the values of all options, including
%     those set by default.
%
% See also: MSAMS.

% Copyright 2010 by Timothy E. Holy

  Y = X;
  options = struct;
  if ~isempty(varargin)
    curarg = 1;
    if isnumeric(varargin{curarg})
      Y = varargin{curarg};
      curarg = curarg+1;
    end
    if (length(varargin) >= curarg)
      options = varargin{curarg};
    end
  end
  options = default(options,'consolidate',true,'n_iter_max',100,'randomize',false);

  %% For each point, compute the closest "uphill" point
  % map0(i) will be the index of the uphill point closest to point i.
  [d,N] = size(X);
  C = eye(d,d);
  M = size(Y,2);
  if (nargout > 1)
    results = struct('n_iter',zeros(1,M),'Ymap',zeros(d,M),'cycling',false(1,M));
  end
  map0 = zeros(1,M);
  n = zeros(1,M);
  for i = 1:M
    isdone = false;
    iter = 0;
    y = Y(:,i);
    ytraj = y;
    ntraj = 0;
    cycling = false;
%     kurt = 0;
    while ~isdone
      if options.randomize
        ri = round(N*rand(1,N)+0.5);
        XX = X(:,ri);
      else
        XX = X;
      end
%       nbrInfo = bn_by_distance(XX,y);
%       y = nbrInfo.nbrhoodMean;
      [ni,nbrLabel] = bn_by_distance_new(XX,y,C,options);
      y = mean(XX(:,nbrLabel(1:ni)),2);
      iter = iter+1;
      [~,minDistIndex] = mindist(y,Y);  % Find the closest probe point
      isdone = (iter >= options.n_iter_max) || (minDistIndex ~= i);
      if ~isdone && ~options.randomize
%         [cycling,ytraj,ntraj,peakIndex] = iscycling(ytraj,ntraj,y,nbrInfo.n);
        [cycling,ytraj,ntraj,peakIndex] = iscycling(ytraj,ntraj,y,ni);
%         kurt(end+1) = nbrInfo.kurtosis; %#ok<AGROW>
      end
      isdone = isdone | cycling;
    end
    map0(i) = minDistIndex;
    if cycling
      n(i) = ntraj(peakIndex);
%       kurt = kurt(peakIndex);
      y = ytraj(:,peakIndex-1);  % Store the one that _produced_ the largest neighborhood
    else
%       n(i) = nbrInfo.n;
       n(i) = ni;
%       kurt = nbrInfo.kurtosis;
    end
    if (nargout > 1)
      results.n_iter(i) = iter;
      results.Ymap(:,i) = y;
      results.cycling(i) = cycling;
%       results.kurtosis(i) = kurt;
    end
  end

  %% Flow everything to the peak
  map = leapfrog(map0,n);

  %% Assign clusters
  [~,~,clust] = unique(map);

  %% Consolidate
  % This is like the code in MSAMS, but here we don't worry about order
  % effects or doing it iteratively (this is a hack)
  if options.consolidate % && isequal(X,Y)
    selfmapIndex = find(map == 1:length(map));
    ops = options;
    ops.randomize = true;
    ops.consolidate = false;
    [clustY,resultsY] = bn_cluster_slow(X,Y(:,selfmapIndex),ops);
    clust = clustY(clust);
    map(selfmapIndex) = selfmapIndex(resultsY.map);
    map = map(map);
  end
%     
%     n_self = length(selfmapIndex);
%     A = zeros(n_self,n_self);
%     for i = 1:length(selfmapIndex)
%       y = Y(:,selfmapIndex(i));
%       nbrInfo = bn_by_distance(X,y);
%       nbrmap = map(nbrInfo.neighborLabel(1:nbrInfo.n));
%       [unbrmap,~,relabelled_nbrmap] = unique(nbrmap);
%       A(i,findainb(unbrmap,selfmapIndex)) = 1;
% %       [~,votes] = agglabel(relabelled_nbrmap);
% %       [~,maxVotesIndex] = max(votes);
% %       map(selfmapIndex(i)) = unbrmap(maxVotesIndex);
%     end
%     comp = connected_components(A);
%     for i = 1:length(comp)
%       thiscomp = selfmapIndex(comp{i});
%       [~,index] = max(n(thiscomp));
%       map(thiscomp) = thiscomp(index);
%     end
%     map = leapfrog(map,n);
%   end
  
  
  if (nargout > 1)
    results.map = map;
    results.map0 = map0;
    results.n = n;
  end