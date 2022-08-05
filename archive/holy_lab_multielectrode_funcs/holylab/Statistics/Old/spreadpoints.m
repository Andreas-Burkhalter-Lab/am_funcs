function [y,ident] = spreadpoints(x,q)
% SPREADPOINTS: choose initial points to occupy all likely clusters
% Syntax:
%   y = spreadpoints(x,q)
%   [y,ident] = spreadpoints(x,q)
% where
%   x is a d-by-N matrix of points in d-dimensional space
%   y is a subset of these points
%   q is the target number of points in the subset
%   ident is an N-vector indexing the closest point in y for each point in
%     x. 
%
% The algorithm starts by drawing random points. Then
%   1. Assign each point in x to its closest representative, defining a
%      set of q "local clusters" (these do not have any real meaning as
%      clusters).
%   2. Compute the width of the longest axis in each "local cluster"
%   3. Iteratively pick the widest "local cluster," and split it by
%      replacing its representative with the two extreme points.
% The idea here is that if there is a real cluster of points which does not
% have a representative, it will add significantly to the width of the
% closest representative(s); by adding the endpoints to the mix you're
% certain to draw a point from the unsampled cluster.
%
% See also: RMEANS.
  
% Copyright 2004 by Timothy E. Holy
  [d,N] = size(x);
  rp = randperm(N);
  y = x(:,rp(1:q)); %  Allocated randomly
  
  % Find the closest y for each x
  [md,ident] = mindist(x,y);
  
  % Find the longest dimension of each "local cluster"
  for i = 1:q
    %[eigval(i),eigvec(:,i)] = largestpc(x(:,find(ident == i)));
    [pc,proj,eigval] = princomp(x(:,find(ident == i))');
    for j = 1:size(proj,2)
      adev(j) = mean(abs(proj(:,j)-median(proj(:,j))));
    end
    lrat{i} = sqrt(eigval')./adev;
  end
  lrat{:}
  lrata = cat(2,lrat{:});
  sort(lrata)
  hist(lrata)
  return
  
  % Iteratively:
  %   Find longest "local cluster" and pick the extreme points of this
  %   cluster as new members of y
  % Quit when the cv of the eigenvalue distribution stops decreasing
  cv = spreadpoints_nzcv(eigval);
  cvnew = cv;
  while (cvnew <= cv(end))
    cv(end+1) = cvnew;
    [mxe,isplit] = max(eigval);
    indx = find(ident == isplit);
    % Calculate projection of points along largest principal component
    xp = eigvec(:,isplit)'*x(:,indx);
    % Replace the "bad" y with end points
    [xps,xpsi] = sort(xp);
    y = [y(:,[1:isplit-1 isplit+1:end]) x(:,indx(xpsi([1 end])))];
    % Recalculate "local cluster" membership and eigenvalue distribution
    [md,ident] = mindist(x,y);
    for i = 1:size(y,2)
      [eigval(i),eigvec(:,i)] = largestpc(x(:,find(ident == i)));
    end
    cvnew = spreadpoints_nzcv(eigval);
  end
 
  function cv = spreadpoints_nzcv(eigval)
    % Coefficient of variation, excluding zeros
    nzindx = find(eigval);
    cv = std(eigval(nzindx))/mean(eigval(nzindx));
    