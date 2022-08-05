function g0 = register_g0(sz,cl)
% REGISTER_G0: create default deformation (identity map)
% This is useful as an initial guess for deformations, when you don't
% have a better guess available.
%
% Syntax:
%   g0 = register_g0(sz)
%   g0 = register_g0(sz,class_string)
% where
%   sz is a 1-by-n_dims vector specifying the size in each dimension;
% and
%   g0 is the "identity map" on that grid, i.e., mapping each point to
%     itself (g(x) = x).  In general, deformations on a grid are
%     specified in the following way:  the grid point defined by
%                        [i,j,....]
%     is mapped to the point
%                        [g{1}(i,j,...), g{2}(i,j,...), ...]
%
% For example:
%   g0 = register_g0([3 4])
% will return a cell array of the following form:
%           [1  1  1  1]
%   g0{1} = [2  2  2  2]
%           [3  3  3  3]
% 
%           [1  2  3  4]
%   g0{2} = [1  2  3  4]
%           [1  2  3  4]
% Consider the grid point [2,3]. We see
%              g0{1}(2,3) = 2, and
%              g0{2}(2,3) = 3,
% so [2,3] is mapped to itself (as are all other grid points).
% A non-trivial deformation would have different values, i.e.,
%           [1.25  1.25  1.25  1.25]
%   g0{1} = [2.25  2.25  2.25  2.25]
%           [3.25  3.25  3.25  3.25]
% 
%           [0.8  1.8  2.8  3.8]
%   g0{2} = [0.8  1.8  2.8  3.8]
%           [0.8  1.8  2.8  3.8]
%           [0.8  1.8  2.8  3.8]
% would encode a rigid translation, where the _coordinates_ move
% "downward" by 0.25 and "leftward" by 0.2.  Note that any function of
% these coordinates (say, an image interpolated using these grid points)
% would therefore "move" in the opposite direction.
  
% Copyright 2006 by Timothy E. Holy
% Help corrected, 2011 by Gary F. Hammen

  if (nargin < 2)
    cl = 'single';
  end
  
  n_dims = length(sz);
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = cast(1:sz(dimIndex),cl);
    x{dimIndex}(end) = cast(sz(dimIndex),cl) - 10*eps(cast(sz(dimIndex),cl)); % prevent NaNs
  end
  g0 = cell(1,n_dims);
  if (n_dims > 1)
    [g0{:}] = ndgrid(x{:});
  else
    g0 = x;
  end
