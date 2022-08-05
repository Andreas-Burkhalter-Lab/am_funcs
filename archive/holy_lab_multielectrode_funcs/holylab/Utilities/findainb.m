function indx = findainb(a,b)
% FINDAINB: find the index of elements of A in B
%
% indx = findainb(a,b)
% a,b may be vectors or cell arrays of strings
% On output, indx is such that a = b(indx).
% Limitation: if b has repeated elements, only the
% index of the first example is used.
%
% See also: INTERSECT, INDEXAINB.

if iscell(b) && ~iscell(a)
    a = {a};
end
% Main challenge: handle repeated elements in a
[sab,index] = sort([a(:);b(:)]);
ia = find(index <= length(a));
% ia contains the indices of the elements in sab that come from a
% Now compute the index of the next element in sab coming from b
% This step is potentially slow if a has a lot of repeated elements
ib = ia+1;
froma = find(index(ib) <= length(a));
while ~isempty(froma)
  ib(froma) = ib(froma)+1;
  %froma = find(index(ib) <= length(a));
  froma = froma(find(index(ib(froma)) <= length(a)));
end
% Check that the values are equal
if ((~iscell(a) & any(sab(ia) ~= sab(ib))) | ...  % This works for vectors
    (iscell(a) & any(~strcmp(sab(ia),sab(ib)))))  % This works for cell array of strings
  error('Not all of a is in b');
end
% Map indices back to original orderings
%isb = index(find(index > length(a))) - length(a);
isb = index(ib) - length(a);
[tmp,iai] = sort(index(ia));
%indx = isb(iai);
indx = isb(iai);
indx = reshape(indx,size(a));
return;

% Old code: didn't handle repeated elements correctly, though
% it certainly is simpler!
[com,ia,ib] = intersect(a(:),b(:));
if (length(com) < length(unique(a)))
  error('Not all of a is in b');
end
[tmp,iai] = sort(ia);   % inverse permutation of a
indx = ib(iai);
