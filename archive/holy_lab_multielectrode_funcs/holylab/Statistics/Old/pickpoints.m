function [xindx,invdb] = pickpoints(x,y,nnew)
% PICKPOINTS: select points in unsampled clusters
% Syntax:
%   xindx = pickpoints(x,y,nnew)
%   [xindx,invdb] = pickpoints(x,y,nnew)
% where
%   x is a d-by-N matrix of points in d dimensions;
%   y is a d-by-q matrix of starter "representative" points;
%   2*nnew is the number of new points you want;
% and
%   xindx is the index of new points, x(:,xindx) being the coordinates.
%
% The algorithm works as follows: each point in x is assigned to its
% closest point in y.  For each such "local cluster" (not a real cluster),
% perform PCA, and compute an index which measures the degree to which
% points are distributed towards the endpoints of each principal component.
% (This measure basically detects dips in the density of
% points projected along a given principal component.) For the nnew local
% clusters & components with the most extreme distribution, return the
% indices of the points at both ends. This should insure that at least one
% of these points resides in an unsampled (or poorly-sampled) cluster.
% 
% The optional return argument invdb is the "inverse dumbbell" value;
% smaller values indicate a greater weighting towards the extremes.
% Standout low values are ones that require splitting; a second call to
% PICKPOINTS can be used to determine if further splitting is needed.
%
% See also: RMEANS.

% Copyright 2004 by Timothy E. Holy

  [d,N] = size(x);
  [d,q] = size(y);
  
  % Find the closest y for each x
  [md,ident] = mindist(x,y);
  
  % For each "local cluster," compute the "dumbbell" index for each
  % principal component (the degree to which points are clustered at the
  % end)
  for i = 1:q
    indx{i} = find(ident == i);
    if (length(indx{i}) > 2*d)
      [pc{i},proj{i},eigval{i}] = princomp(x(:,indx{i})');
      for j = 1:d
        adev(j) = mean(abs(proj{i}(:,j)-median(proj{i}(:,j))));
      end
      % Zero-pad eigenvalues for those cases where the number of points in
      % a local cluster is less than d, if this ever happens.
      % (It may not depending on the value chosen for the number of points
      % in the local cluster.)
      eigval{i} = [eigval{i}' zeros(1,d-length(eigval{i}))];
      dbindx{i} = adev./sqrt(eigval{i});  % Small values are more dumbbell-like
      % (yes, this seems inverted, but it handles NaNs gracefully this way)
    else
      dbindx{i} = nan(1,d);
    end
  end
  
  % For the chosen number of directions, pick the points at opposite ends
  % of the dumbbell and return their indices
  dbindxa = cat(2,dbindx{:});
  [sdb,sdbi] = sort(dbindxa);
  if any(isnan(sdb(1:nnew)))
    warning('Some points will be randomly chosen')
  end
  sdbi = sdbi(1:nnew);           % Just use the most-dumbbell-like ones
  cnum = ceil(sdbi/d);           % "Local cluster" number
  pcnum = sdbi - (cnum-1)*d;     % principal component number within loc. cl.
  xindx = zeros(1,2*nnew);
  for i = 1:nnew
    [minp,minpi] = min(proj{cnum(i)}(:,pcnum(i)));
    [maxp,maxpi] = max(proj{cnum(i)}(:,pcnum(i)));
    xindx(2*i-1:2*i) = indx{cnum(i)}([minpi maxpi]);
  end
  if nargout > 1
    invdb = sdb;
  end
