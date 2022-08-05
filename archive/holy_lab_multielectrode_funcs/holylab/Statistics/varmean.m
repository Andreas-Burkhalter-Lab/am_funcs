function [v,m] = varmean(X,flag,dim)
% VARMEAN: Compute variance and mean & save memory
% Syntax:
%   [v,m] = varmean(X,flag,dim)
% where
%   X is an array;
%   flag if flag is true, normalizes by n, otherwise by n-1
%   dim is the dimension along which the variance is computed (default:
%      first non-singleton dimension)
% and
%   v is the variance array;
%   m is the mean array.
%
% See also: MEAN, VAR.
   
% Copyright 2004 by Timothy E. Holy

% Stuff cut from old version:
% This routine loops over the coordinate in the chosen dimension; this is
% slower than the algorithm used in VAR unless you find yourself short on
% memory. Plus, it returns both the mean and the variance.
%
%    for i = 1:ndims(X)
%       indxstr{i} = ':';
%    end
%    v = zeros(size(m));
%    for i = 1:n
%       indxstr{dim} = i;
%       xshift = double(reshape(X(indxstr{:}),size(m))) - m;
%       v = v+xshift.^2;
%    end

   if (nargin < 2)
     flag = 0;
   end
   if (nargin < 3)
      dim = find(size(X) ~= 1,1);
      if isempty(dim)
         dim = 1;
      end
   end
   n = size(X,dim);
   if (n < 2 && flag == 0)
      error('Variance requires more than one replicate');
   end
   m = sum(X,dim)/n;
   dX = bsxfun(@minus,X,m);
   v = sum(dX.^2,dim);
   if flag
     v = v/n;
   else
     v = v/(n-1);
   end
