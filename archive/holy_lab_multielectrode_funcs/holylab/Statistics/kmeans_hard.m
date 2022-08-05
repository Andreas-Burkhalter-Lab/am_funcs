function [idx,c,rmsd] = kmeans_hard(X,k,options)
% KMEANS_HARD: fast "hard K-means" algorithm
% Syntax:
%   [idx,c,rmsd] = kmeans_hard(X,k)
%   [idx,c,rmsd] = kmeans_hard(X,c0)
%   ... = kmeans_hard(...,options)
% where
%   X is your d-by-npts matrix of data points (note this is transpose wrt
%     Matlab's data matrix)
%   k is the number of clusters you want
%                OR
%   c0 is the starting guess for the cluster centers, a matrix of size
%     d-by-nclusts
%   options is a structure which may have the following fields:
%     showprogress (default false): if true, prints the number of points
%       changing assignment on each iteration.
%     remove_empty (default true): if true, empty clusters are not returned
% and
%   idx is the cluster number for each point (1-by-npts)
%   c is a d-by-nclusts matrix giving the cluster centers
%   rmsd is a 1-by-nclusts vector with the root-mean-square distance of
%     points from the cluster center
%
% This algorithm is much faster than Matlab's kmeans, but only does the
% "batch" updates, and skips the "online" updates. (However, it's faster
% even on the batch updates, largely because of using the MEX file
% mindist.)
%
% See also: EM_GAUSS, KMEANS.
  
% Copyright 2005 by Timothy E. Holy
  
  [d,npts] = size(X);
  if isscalar(k)
    startindx = randperm(npts);
    c = X(:,startindx(1:k));
  else
    c = k;
    k = size(c,2);
    if (size(c,1) ~= d)
      error('Dimensionality mismatch between c0 (initial cluster centers) and data');
    end
  end
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'remove_empty',true,'showprogress',false,'plotprogress',false,'pauseeach',false);
  
  if (options.plotprogress && d >= 2)
    options = kmeans_hard_plot(X,c,options);
  end

  idxold = zeros(1,npts);  % Cluster assignment on previous iteration
  [md,idx] = mindist(X,c); % Cluster assignment & distance on current iter.
  moved = true;     % Flag to see if any points have changed assignment
  iter = 0;
  while moved
    clabel = agglabel(idx);
    if (length(clabel) < k)
      clabel{k} = [];
    end
    moved = false;
    for i = 1:k
      if any(idx(clabel{i}) ~= idxold(clabel{i}))
        moved = true;
        c(:,i) = mean(X(:,clabel{i}),2);  % Recompute cluster center position
      end
    end
    if options.showprogress
      fprintf('  %d\n',sum(idx ~= idxold))
    end
    if (options.plotprogress && d >= 2)
      options = kmeans_hard_plot(X,c,options);
    end
    idxold = idx;
    [md,idx] = mindist(X,c);  % Reassign points
    iter = iter+1;
  end
  %iter
  if options.remove_empty
    l = cellfun(@length,clabel);
    keepFlag = l > 0;
    clabel = clabel(keepFlag);
    c = c(:,keepFlag);
    cskf = cumsum(keepFlag);
    idx = cskf(idx);
  end
  % If requested, supply rms sizes
  if (nargout > 2)
    rmsd = zeros(1,k);
    for i = 1:k
      rmsd(i) = sqrt(mean(md(clabel{i})));
    end
  end
  
function options = kmeans_hard_plot(X,c,options)
  plot(X(1,:),X(2,:),'.');
  line(c(1,:),c(2,:),'LineStyle','none','Marker','x','Color','r','MarkerSize',16,'LineWidth',2);
  axis equal
  if options.pauseeach
    pause
    c = get(gcf,'CurrentCharacter');
    if (c == 'q')
      options.pauseeach = false;
    end
  else
    drawnow
  end
    