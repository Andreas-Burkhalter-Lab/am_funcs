function bricks = prepare_bricks(arraysz,bricksz,brickoffset)
% PREPARE_BRICKS: divide array (or image) volume into bricks of common size
% Syntax:
%   bricks = prepare_bricks(arraysz,bricksz)
%   bricks = prepare_bricks(arraysz,bricksz,brickoffset)
% where all of the inputs are size vectors.  The first two contain the
% size of the  array and the size of the brick, respectively. The last is
% an optional offset, if you don't want the first brick to start at
% 1,1,...
% The output is a cell array, each element of which is a cell array
% containing the coordinates for the corresponding brick.
%
% Example
%   bricks = prepare_bricks(size(im),[5 5]);
% would prepare a list of 5x5 bricks.  The ith brick can be extracted as
%   x = brick{i};
%   imsnip = im(x{:});

% Copyright 2010 by Timothy E. Holy
  
  n_dims = length(arraysz);
  if (nargin < 3)
    brickoffset = zeros(1,n_dims);
  end
  xstart = cell(1,n_dims);
  xend = cell(1,n_dims);
  for dimIndex = 1:n_dims
    xstart{dimIndex} = brickoffset(dimIndex)+1:bricksz(dimIndex): ...
      max(arraysz(dimIndex)-bricksz(dimIndex)/4,brickoffset(dimIndex)+1);
    xend{dimIndex} = [xstart{dimIndex}(2:end)-1 arraysz(dimIndex)];
  end
  Xstart = cell(1,n_dims);
  Xend = cell(1,n_dims);
  [Xstart{:}] = ndgrid(xstart{:});
  [Xend{:}] = ndgrid(xend{:});
  n_bricks = numel(Xstart{1});
  bricks = cell(1,n_bricks);
  x = cell(1,n_dims);
  for brickIndex = 1:n_bricks
    for dimIndex = 1:n_dims
      x{dimIndex} = Xstart{dimIndex}(brickIndex): ...
        Xend{dimIndex}(brickIndex);
    end
    bricks{brickIndex} = x;
  end
end
