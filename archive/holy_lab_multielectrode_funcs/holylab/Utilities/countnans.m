function s = countnans(A,dim)
% countnans: count the number of NaNs in coordinate slices
%
% Syntax:
%   s = countnans(A,dim)
% where
%   A is a multidimensional array
%   dim is the dimension along which you want to slice A
% and
%   s is a vector, s(i) is the total number of NaNs in the slice
%     A(:,...,:,i,:,...,:) (where the i is in the "dim"th coordinate).
%
% Example: let
%    A = [nan(1,4) 4; nan 1 2 3 4; nan nan 2 3 4; nan 1 2 3 nan]
% Then
%    countnans(A,1)
% returns [4; 1; 2; 2]
% and
%    countnans(A,2)
% returns [4 2 1 1 1].

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    dim = 1;
  end
  s = isnan(A);
  for j = ndims(A):-1:1
    if (j ~= dim)
      s = sum(s,j);
    end
  end
