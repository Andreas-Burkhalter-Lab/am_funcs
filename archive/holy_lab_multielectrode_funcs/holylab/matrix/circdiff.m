function dA = circdiff(A,dimIndex,mode)
% circdiff: finite differences with periodic boundary conditions
%
% This function produces finite differences (an estimate of the derivative)
% under the assumption of periodic boundary conditions. Thus, unlike the
% ordinary "diff", the output arrays have the same size as the input array.
%
% Syntax:
%   dA = circdiff(A,dimIndex)
%   dA = circdiff(A,dimIndex,mode)
% where
%   A is an array of arbitrary dimensionality
%   dimIndex is the dimension along which you want to compute the finite
%     difference
%   mode (optional) is a string which controls the sense of the derivative:
%     in one dimension, with 'forward' (the default) one has 
%        dA(i) = A(i+1)-A(i)
%     and with 'reverse' or 'backward' one has
%        dA(i) = A(i) - A(i-1)
%
% See also: diff.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 3)
    mode = 'forward';
  end
  x = repmat({':'},1,ndims(A));
  xend = x; xend{dimIndex} = size(A,dimIndex);
  x1 = x; x1{dimIndex} = 1;
  dA1 = A(x1{:}) - A(xend{:});
  switch mode
    case 'forward'
      dA = diff(A,1,dimIndex);
      dA(xend{:}) = dA1;
    case {'reverse','backward'}
      dA = zeros(size(A),class(A));
      xa = x; xa{dimIndex} = 2:size(A,dimIndex);
      dA(xa{:}) = diff(A,1,dimIndex);
      dA(x1{:}) = dA1;
    otherwise
      error('mode not recognized');
  end
end
