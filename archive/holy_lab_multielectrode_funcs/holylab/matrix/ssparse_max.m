function s = ssparse_max(M,dimIndex)
% SSPARSE_MAX: compute maximum over rows or columns of "sparse structure matrices"
% Syntax:
%   s = ssparse_max(M,dimIndex)
% where
%   M is a sparse structure matrix
%   dimIndex is 1 or 2
% and
%   s is a row or column vector yielding the max of M over the given
%     dimension.

% Copyright 2010 by Timothy E. Holy

  sz = size(M);
  if (sz(dimIndex) == 1)
    % This is the easy case, we are taking the max along the sorting
    % dimension
    s = arrayfun(@(p) max([p.value 0]),M);
  else
    if (sz(1) == 1)
      % Column-major ordering, taking max over columns
      index = vertcat(M.rowi);
      value = vertcat(M.value);
      s = maxlabel(index,value);
    else
      % Row-major ordering, summing over rows
      index = horzcat(M.coli);
      value = horzcat(M.value);
      s = maxlabel(index(:),value(:))';
    end
  end
end
      
function s = maxlabel(index,value)
  clabel = agglabel(index);
  n_labels = length(clabel);
  s = -inf(n_labels,1);
  for i = 1:n_labels
    if ~isempty(clabel{i})
      s(i) = max(value(clabel{i}));
    end
  end
end