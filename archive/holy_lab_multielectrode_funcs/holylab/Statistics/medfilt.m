function y = medfilt(x,varargin)
% MEDFILT: median filtering with proper end treatment
% Syntax: see medfilt1
%
% medfilt1 pads the ends with zeros; this algorithm then fixes up the ends.
% The points contaminated by the padding are replaced by the value obtained
% from the largest-sized median filter which does not need padding.
% Therefore, an end point is copied unmodified, the next point in is the
% median of the 3 points at the edge, and so forth.
%
% Strongly recommended:  use only odd-sized filters.
%
% Todo: support multidimensional data
%
% See also: MEDFILT1.

% Copyright 2004 by Timothy E. Holy
y = medfilt1(x,varargin{:});
if (min(size(x)) > 1)
  error('multidimensional data not yet supported');
end
if ~isempty(varargin)
  n = varargin{1};
else
  n = 3;
end
rem = mod(n+1,2);
padsize = min(ceil((n-1)/2),ceil(length(x)/2));
for i = 1:padsize
  rngsize = 2*i-1+rem;
  y(i) = median(x(1:rngsize));
  y(end-i+1) = median(x(end-rngsize+1:end));
end