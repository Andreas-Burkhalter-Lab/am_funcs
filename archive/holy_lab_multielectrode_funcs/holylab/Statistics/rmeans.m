function [clust,nclust] = rmeans(x,options)
% RMEANS: Cluster points non-parametrically
% This works by a modified mean-shift algorithm, in that it
%   (1) Picks a set of landmark points (by kmeans);
%   (2) Assigns each input point to the closest landmark;
%   (3) Moves each landmark "uphill," i.e., up the density gradient of
%       the input points, and merges the ones that end up at the same
%       position into a single cluster
%
% Syntax:
%   clust = rmeans(x)
%   clust = rmeans(x,options)
% where
%   x is a d-by-N matrix of points in d-dimensional space;
%   options is an options structure with the following fields:
%     plotclust (optional): if true, will show the clusters (in the first
%       two dimensions) at the end of clustering
%     plotproject (optional): if supplied, will do all visualization in
%       terms of x coordinates projected on these vectors. Should be a
%       2-by-d matrix.
%     q (default 30): number of landmarks to use
%     Rfac (default 2): a factor used to adjust the initially-chosen R
%       for less (Rfac < 1) or more (Rfac > 1) smoothing.
%       R = Rfac*Rinitial (Rinitial is the root-mean-square
%       distance of points from each landmark). 
%  + any options valid for MEANSHIFT.
%
%  On output,
%    clust is a 1-by-N vector, giving the cluster number assigned to each
%      point in x
%
% See also: MEANSHIFT.
  
% Copyright 2004 Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  options = rmeansoptions(options);
  [d,N] = size(x);
  if ~isfield(options,'plotproject')
    options.plotproject = zeros(2,d);
    options.plotproject([1 2],[1 2]) = eye(2);
  end

  % Choose landmarks & assign points to closest landmark
  [yclose,y,R] = kmeans_hard(x,options.q);
  % Adjust R if needed
  if (0)  % The following code is never executed, but can be useful for
          % debugging/inspection, so I left it in place
    % Distance to nearest landmark (other than self)
    ny = size(y,2);
    for i = 1:ny;
      iother = setdiff(1:ny,i);
      dy = repmat(y(:,i),1,ny-1)-y(:,iother);
      Rd(i) = sqrt(min(sum(dy.^2,1)));
    end
  end
  nzindx = R > 0;
  R(~nzindx) = mean(R(nzindx)); % set any with R=0 to mean of nonzero
  R = R*options.Rfac;

  % Do the meanshift algorithm
  clust = meanshift(x,y,R,options);
  % Map cluster numbers from landmarks to points  
  clust = clust(yclose);
  
  % Plot, if desired
  if options.plotclust
    clf
    co = get(gca,'ColorOrder');
    hold on
    clabel = agglabel(clust);
    nclust = length(clabel);
    for i = 1:nclust
      colindx = mod(i-1,size(co,1))+1;
      clindx = clabel{i};
      xp = options.plotproject*x;
      if (size(xp,1) > 1)
        plot(xp(1,clindx),xp(2,clindx),'.','Color',co(colindx,:));
      else
        plot(xp(clindx),zeros(size(clindx)),'.','Color',co(colindx,:));
      end
    end
    hold off
    axis equal
  end
  
  
function op = rmeansoptions(op)
  if ~isfield(op,'ploteach')
    op.ploteach = 0;
  end
  if ~isfield(op,'plotclust')
    op.plotclust = 0;
  end
  if ~isfield(op,'shrinkfrac')
    op.shrinkfrac = 0.3;
  end
  if ~isfield(op,'q')
    op.q = 30;
  end
  if ~isfield(op,'Rfac')
    op.Rfac = 2;  % This seems to work surprisingly well across many
                  % conditions!
  end