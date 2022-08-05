function indx = time2indx(t,tref,varargin)
% TIME2INDX: find the index of an event occuring near particular times
% Syntax:
%   indx = time2indx(t,tref)
%   indx = time2indx(t,tref,dt)
%   indx = time2indx(...,'sorted')
% where
%   t is a vector of approximate times at which particular events
%     occured;
%   tref is the vector of exact times at which all events occured;
%   dt (optional) is the maximum tolerance you'll allow in matching
%     approximate and exact times;
% and
%   indx is a vector of integers such that tref(indx) is the set of exact
%     times closest to the user-specified times in t.  If no match can be
%     found within a time dt, the corresponding indx value is NaN.
%
% If the string 'sorted' is appended to the argument list, the algorithm
% will assume that the input data are sorted; this speeds it greatly.

% Copyright 2004 Timothy E. Holy
  
flag = 'unsorted';
if (nargin < 3)
  dt = Inf;
else
  if (length(varargin) > 1 | isnumeric(varargin{1}))
    dt = varargin{1};
  end
  if (ischar(varargin{end}))
    flag = varargin{end};
  end
end
% Sort the times to make it efficient
%( O(length(tref)) + O(length(t)) rather than their product)
if ~strcmp(flag,'sorted')
  [stref,strefi] = sort(tref);
  [st,sti] = sort(t);
else
  stref = tref;
  st = t;
end
nt = length(t);
ntref = length(tref);
indx = nan(1,nt);
j = 1;
for i = 1:nt
  dtbest = abs(stref(j)-st(i));
  while (j < ntref && abs(stref(j+1)-st(i)) < dtbest)
    dtbest = abs(stref(j+1)-st(i));
    j = j+1;
  end
  if (dtbest < dt)
    indx(i) = j;
  end
end
% Now map the index back to the original ordering
if ~strcmp(flag,'sorted')
  % First compensate for the sorting of t
  [ssti,stiinv] = sort(sti);
  indx = indx(stiinv);
  % Now compensate for the sorting of tref
  nnan = find(~isnan(indx));
  indxo = nan(size(indx));
  indxo(nnan) = strefi(indx(nnan));
  indx = indxo;
end