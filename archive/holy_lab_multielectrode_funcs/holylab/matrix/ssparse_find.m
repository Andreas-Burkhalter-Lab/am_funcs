function [index1,index2,value] = ssparse_find(M)
% SSPARSE_FIND: extract indices and values of a "sparse structure matrix"
% Syntax:
%   [index1,index2,value] = ssparse_find(M)
% The sparse matrix S corresponding to M would be specified as
%     S = sparse(index1,index2,value)

% Copyright 2009 by Timothy E. Holy

  if (size(M,1) == 1)
    % Column major ordering
    index1 = vertcat(M.rowi);
    index2 = make_counting_index(arrayfun(@(p) length(p.rowi),M))';
    value = vertcat(M.value);
  else
    % Row major ordering
    index1 = make_counting_index(arrayfun(@(p) length(p.coli),M))';
    index2 = horzcat(M.coli);
    value = horzcat(M.value);
  end
end