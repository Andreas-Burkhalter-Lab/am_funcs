function dA = deriv(A,dimIndex,order)
% DERIV: estimate the derivative of a matrix along a coordinate
% Unlike DIFF, this returns a matrix of the same size.
% Syntax:
%   dA = deriv(A,dimIndex)
%   dA = deriv(A,dimIndex,order)
% where
%   A is an array (can be multidimensional)
%   dimIndex is the dimension along which the derivative should be taken
%   order (default 1) is the order of the derivative (must be 1 or 2)
% and
%   dA is the discrete derivative of A along coordinate dimIndex.
% 
% This function uses centered differencing to get a higher-order estimate
% of the derivative.
% This is very much like GRADIENT, except:
%   (1) It computes the derivative only along the requested coordinate;
%   (2) It doesn't swap the first two dimensions. Thus it may be better
%       suited for multidimensional problems.
% It's also slightly faster, since it doesn't call PERMUTE.
%
% See also: GRADIENT, DIFF.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    order = 1;
  end
  if (order ~= 1 && order ~= 2)
    error('order must be 1 or 2');
  end
  clsA = class(A);
  if (clsA(1) == 'u')
    error('deriv does not support unsigned types');
  end
  
  sz = size(A);
  if (sz(dimIndex) == 1)
    error('Can''t compute the derivative in a dimension with size 1');
  end
  n_dims = ndims(A);
  colons = repmat({':'},1,n_dims);
  dA = zeros(sz,class(A));
  
  colons_p = colons;
  colons_m = colons;
  colons_c = colons;
  if (order == 1)
    % For the middle of the matrix, do centered differencing
    if (sz(dimIndex) > 2)
      colons_p{dimIndex} = 3:sz(dimIndex);
      colons_m{dimIndex} = 1:sz(dimIndex)-2;
      colons_c{dimIndex} = 2:sz(dimIndex)-1;
      dA(colons_c{:}) = (A(colons_p{:}) - A(colons_m{:}))/2;
    end
    % For the edges, use the one-sided derivatives
    colons_c{dimIndex} = 1;
    colons_p{dimIndex} = 2;
    dA(colons_c{:}) = A(colons_p{:}) - A(colons_c{:});
    colons_c{dimIndex} = sz(dimIndex);
    colons_m{dimIndex} = sz(dimIndex)-1;
    dA(colons_c{:}) = A(colons_c{:}) - A(colons_m{:});
  else
    error('2nd order derivative not yet implemented');
  end
  
  