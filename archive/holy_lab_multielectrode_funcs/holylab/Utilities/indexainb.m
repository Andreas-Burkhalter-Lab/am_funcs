function indx = indexainb(a,b,sorted)
% indexainb: find the index of those elements of A that appear in B
%
% indx = indexainb(a,b)
% a,b may be vectors or cell arrays of strings
% On output, indx is such that a(indx) are elements of b
%
% indx = indexainb(a,b,'sorted') assumes that
% a and b are already in monotonically increasing order,
% and is faster in that case. However, it will give
% you incorrect answers if you are wrong about this.
%
% See also: INTERSECT, FINDAINB.

% Copyright Timothy E. Holy <holy@pcg.wustl.edu>

% The workhorse is sort (very fast), though we call
% it up to 4 times!
% Sort preserves the order of repeated elements.
% Therefore, by sorting in a;b and then b;a, we
% can find out which elements didn't go past each other
presort = 1;
if (nargin > 2 & ischar(sorted) & strcmp(sorted,'sorted'))
    [sab,index1] = sort([a(:);b(:)]);
    [sab,index2] = sort([b(:);a(:)]);
    presort = 0;
else
    [sa,aindex] = sort(a(:));
    [sb,bindex] = sort(b(:));
    [sab,index1] = sort([sa;sb]);
    [sab,index2] = sort([sb;sa]);
end
% Focus on the a elements
ia = find(index1 <= length(a));
% Convert the indices in the b;a sort to
% indices as if it had been an a;b sort
newib = index2(ia);
iib1 = find(newib <= length(b));
iib2 = find(newib > length(b));
newib(iib1) = newib(iib1) + length(a);
newib(iib2) = newib(iib2) - length(b);
% Now, any index which is different between sorts is
% a keeper
indx = find(index1(ia) ~= newib);
if presort
    indx = sort(aindex(indx));
end
