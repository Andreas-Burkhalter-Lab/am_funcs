function L = laplacian(sz)
% LAPLACIAN: calculate the laplacian operator on a grid
% This laplacian operator is defined so that, for two functions f and g
% defined on the grid,
%    sum (gradf dot gradg) = - (f Lg)
% where Lg is the Laplacian of g.
%
% The ordering of grid points in L is the native Matlab ordering for
% arrays, i.e., first columns, then rows, then ...
%
% Syntax:
%   L = laplacian(sz)
% where
%   sz is a vector giving the number of grid points along each dimension;
% and
%   L is a sparse matrix defining the Laplacian operator.
  
  n_dims = length(sz);
  colons = repmat({':'},1,n_dims);
  n_pts = prod(sz);
  pt_index = reshape(1:n_pts,sz);
  L = sparse([],[],[],n_pts,n_pts,(2*n_dims+1)*n_pts);
  for dimIndex = 1:n_dims
    bwIndex = colons;
    fwIndex = colons;
    bwIndex{dimIndex} = 1:sz(dimIndex)-1;
    fwIndex{dimIndex} = 2:sz(dimIndex);
    bwPair = pt_index(bwIndex{:});
    fwPair = pt_index(fwIndex{:});
    indexL = sub2ind(size(L),bwPair(:),fwPair(:));
    L(indexL) = 1;
  end
  L = L + L';
  sumL = full(sum(L));
  L = L - spdiags(sumL',0,n_pts,n_pts);
  
  