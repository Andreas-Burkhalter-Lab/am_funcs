function npb = HistSplit(t,tsplit)
% HistSplit: histogram, given the bin boundaries
% Syntax:
%   npb = HistSplit(t,tsplit)
% where
%   t contains the data points;
%   tsplit contains the bin boundaries;
% and npb is a vector of size length(tsplit)+1.
%   npb(1) is the number of data points less than tsplit(1);
%   npb(2) is the number of data points between tsplit(1) and tsplit(2);
%   ...
%   npb(end) is the number of data points larger than tsplit(end).
%
% See also: SPLITEVENLY.
  
% Copyright 2001 by Timothy E. Holy
  
nsplitm1 = length(tsplit) - 1;
[ts,ii] = sort([t(:)',tsplit(:)']);
indx = find(ii > length(t));
cnpb = indx(end-nsplitm1:end) - (0:nsplitm1);
npb = diff([1,cnpb,length(t)+1]);
