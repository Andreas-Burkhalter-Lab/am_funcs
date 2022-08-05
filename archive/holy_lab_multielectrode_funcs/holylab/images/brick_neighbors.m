function varargout = brick_neighbors(bricks)
% BRICK_NEIGHBORS: find the bricks that share vertices
%
% Syntax:
%   isnbr = brick_neighbors(bricks)
%   [nbr1,nbr2] = brick_neighbors(bricks)
% where
%   bricks is a cell array of the type output by prepare_bricks;
%   bricksz is a 1-by-n_dims vector giving the number of pixels along each
%     axis of each brick
% and
%   isnbr is a n_bricks-by-n_bricks matrix, where isnbr(i,j) is true if i
%     and j are neighbors
%  OR
%   nbr1 and nbr2 are vectors, where nbr1(i) is a neighbor of nbr2(i).
%     Each pair is listed only once, and no brick is a neighbor of itself.
%
% See also: PREPARE_BRICKS.

% Copyright 2010 by Timothy E. Holy

  n_bricks = length(bricks);
  n_dims = length(bricks{1});
  % Collect "upper-left" vertex and infer the bricksz
  X = zeros(n_dims,n_bricks);
  sz = zeros(n_dims,n_bricks);
  for i = 1:n_bricks
    for j = 1:n_dims
      X(j,i) = bricks{i}{j}(1);
      sz(j,i) = length(bricks{i}{j});
    end
  end
  bricksz = max(sz,[],2);
  Xnorm = X ./ repmat(bricksz(:),1,n_bricks);
  isnbr = true(n_bricks,n_bricks);
  for j = 1:n_dims
    sd = sqrdist(Xnorm(j,:),Xnorm(j,:));
    isnbr = isnbr & (sd < 1.01);
  end
  if (nargout == 1)
    varargout = {isnbr};
  else
    [nbr1,nbr2] = find(isnbr);
    keepFlag = nbr1 > nbr2;
    nbr1 = nbr1(keepFlag);
    nbr2 = nbr2(keepFlag);
    varargout = {nbr1,nbr2};
  end
  