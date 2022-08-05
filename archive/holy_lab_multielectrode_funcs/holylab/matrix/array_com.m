function [x,xvar] = array_com(A)
% array_com: compute center of mass of an array
% Syntax:
%   x = array_com(A)
% where
%   A is an array of arbitrary dimensionality
% and
%   x is a 1-by-n_dims vector containing the center of mass along each
%     coordinate.
%
%   [x,xvar] = array_com(A)
% also computes the coordinate variance of the array.
%
% Note that if A has negative entries, the center of mass is an ill-defined
% concept. The non-negativity of A is not enforced by this function, so if
% you want non-negative A you need to arrange this yourself.

% Copyright 2010 by Timothy E. Holy

  sz = size(A);
  n_dims = ndims(A);
  cl = class(A);
  x = zeros(1,n_dims,cl);
  if (nargout > 1)
    xvar = zeros(1,n_dims,cl);
  end
  Asum = sum(A(:));
  o = ones(1,n_dims,cl);
  for dimIndex = 1:n_dims
    dimvec = o;
    dimvec(dimIndex) = sz(dimIndex);
    Ax = bsxfun(@times,A,reshape(1:sz(dimIndex),dimvec));
    x(dimIndex) = sum(Ax(:))/Asum;
    if (nargout > 1)
      Ax = bsxfun(@times,A,reshape(((1:sz(dimIndex))-x(dimIndex)).^2,dimvec));
      xvar(dimIndex) = sum(Ax(:))/Asum;
    end
  end
end
