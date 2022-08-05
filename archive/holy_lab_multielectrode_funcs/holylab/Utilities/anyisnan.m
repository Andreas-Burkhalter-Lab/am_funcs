function flag = anyisnan(A,dim)
% anyisnan: check whether any entry in coordinate slices is NaN
%
% Syntax:
%   flag = anyisnan(A,dim)
% where
%   A is a multidimensional array
%   dim is the dimension along which you want to slice A
% and
%   flag indicates whether any entries are NaN in each slice of A; flag(i)
%     is true if any of A(:,...,:,i,:,...,:) is NaN (where the i is in the
%     "dim"th coordinate).
%
% Example: let
%    A = [1 2 3; NaN 5 6];
% Then
%    anyisnan(A,1)
% returns [false; true]
% and
%    anyisnan(A,2)
% returns [true false false].

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    dim = 1;
  end
  flag = isnan(A);
  for j = ndims(A):-1:1
    if (j ~= dim)
      flag = any(flag,j);
    end
  end
