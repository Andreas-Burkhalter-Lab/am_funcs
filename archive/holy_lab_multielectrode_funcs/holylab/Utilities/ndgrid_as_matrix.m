function gridout = ndgrid_as_matrix(varargin)
% NDGRID_AS_MATRIX: give grid coordinates in column format
% Syntax:
%   gridout = ndgrid_as_matrix(x1,x2,x3,...)
% where
%   x1,x2,x3,... are vectors, usually of coordinate values;
% and
%   gridout's rows contain all possible combinations of the input vector
%     values.
%
% See also: NDGRID.
  
% Copyright 2006 by Timothy E. Holy
  
  n_dims = length(varargin);
  if (n_dims == 1)
    gridout = varargin{1}(:);
    return
  end
  gridtmp = cell(1,n_dims);
  [gridtmp{:}] = ndgrid(varargin{:});
  for i = 1:n_dims
    gridtmp{i} = gridtmp{i}(:);
  end
  gridout = cat(2,gridtmp{:});
  