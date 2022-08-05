function index = make_counting_index(l)
% MAKE_COUNTING_INDEX: create a vector to index variable-length blocks
%
% When concatenating structure or cell arrays, you often want an index back
% to the original object. For example, suppose C is a cell array of
% vectors, and l(i) is the length of C{i}. (You can create l with l =
% cellfun(@length,C).) Then when you do
%    c = cat(2,C{:}),
% you can determine the element of C corresponding to each element in c in
% this way:
%    index = make_counting_index(l).
% In other words, c(k) came from the index(k)-th element of C.

% Copyright 2009 by Timothy E. Holy

  lcum = cumsum(l(:));
  index = accumarray(lcum(:)+1,ones(size(lcum(:))));
  index = cumsum(index)+1;
  index = index(1:end-1);
  if (size(l,1) == 1)
    index = index';
  end
end
  
  
