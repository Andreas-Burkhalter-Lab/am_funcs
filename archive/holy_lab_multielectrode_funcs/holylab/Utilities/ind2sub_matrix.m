function sub = ind2sub_matrix(sz,ind,test)
% ind2sub_matrix: like ind2sub, but coordinates are specified as a matrix
% Syntax:
%   sub = ind2sub_matrix(sz,ind)
% where
%   sz is a vector holding the size of the array
%   ind is the linear index into the array, a vector of length n_pts
% and
%   sub is a n_pts-by-n_dims matrix, where each row holds the coordinates
%     of the corresponding point in the array
%
%   sub = ind2sub_matrix(sz,ind,test)
% allows you to control whether the index is tested for being in-bounds. If
% test is false, the test is omitted. This saves time. Default is
% test=true.
%
% See also: ind2sub, sub2ind_matrix.
  
% Copyright 2011 by Timothy E. Holy

  index_skip = cumprod([1 sz]);
  if (nargin < 3 || test)
    if any(ind < 1) || any(ind > index_skip(end))
      error('Index is out of range');
    end
  end
  n_dims = length(sz);
  sub = zeros(length(ind),n_dims);
  ind = ind(:)-1;
  for dimIndex = n_dims:-1:1
    x = floor(ind/index_skip(dimIndex));
    sub(:,dimIndex) = x+1;
    ind = ind-x*index_skip(dimIndex);
  end
  