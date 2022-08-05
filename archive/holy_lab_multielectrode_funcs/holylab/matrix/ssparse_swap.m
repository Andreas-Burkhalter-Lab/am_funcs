function [Ms,index1,index2,value] = ssparse_swap(M,max_missing)
% SSPARSE_SWAP: change from column-major to row-major ordering (or vice versa)
%
% This function takes a column-major "sparse structure array" M, where M is
% of size 1-by-n and has the row index field (rowi), and returns a sparse
% structure array Ms of size m-by-1 which has a column index field (coli).
% TODO: It also implements the reverse transformation.
%
% Syntax:
%   Ms = ssparse_swap(M)
%   Ms = ssparse_swap(M,m_min)
% This variant will insure that the output array has a least m_min elements
% (conceptually, m_min is the minimum number of rows in the output)

%   [Ms,index1,index2,value] = ssparse_swap(...)
% This variant returns vectors that would be equivalent to
%     [index1,index2,value] = find(S)
% for a sparse matrix S.

% Copyright 2009 by Timothy E. Holy

  [index1,index2,value] = ssparse_find(M);
  if (size(M,1) == 1)
    % Column major ordering, switch to row ordering
    [index1s,sortOrder] = sort(index1);
    breaks = [0;find(diff(index1s) > 0);length(index1)];
    n = index1s(end);
    if (nargin > 1)
      n = max(n,max_missing);
    end
    Ms = repmat(struct('coli',[],'value',[],'coli_min',1,'coli_max',length(M)),n,1);
    index2 = index2';
    value = value';
    for k = 1:length(breaks)-1
      rng = breaks(k)+1:breaks(k+1);
      sortOrderTmp = sortOrder(rng);
      thisindx = index1s(rng(1));
      Ms(thisindx).coli = index2(sortOrderTmp);
      Ms(thisindx).value = value(sortOrderTmp);
    end
  end
