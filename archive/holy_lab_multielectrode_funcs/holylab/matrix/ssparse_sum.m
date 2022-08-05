function s = ssparse_sum(M,dimIndex)
% SSPARSE_SUM: compute sums over rows or columns of "sparse structure matrices"
% Syntax:
%   s = ssparse_sum(M,dimIndex)
% where
%   M is a sparse structure matrix
%   dimIndex is 1 or 2
% and
%   s is a row or column vector yielding the sum of M over the given
%     dimension.

% Copyright 2009 by Timothy E. Holy

  sz = size(M);
  if (sz(dimIndex) == 1)
    s = arrayfun(@(p) sum(p.value),M);
  else
    if (sz(1) == 1)
      % Column-major ordering, summing over columns
      index = vertcat(M.rowi);
      value = vertcat(M.value);
      s = accumarray(index,value);
    else
      % Row-major ordering, summing over rows
      index = horzcat(M.coli);
      value = horzcat(M.value);
      s = accumarray(index(:),value(:))';
    end
  end
end
      
