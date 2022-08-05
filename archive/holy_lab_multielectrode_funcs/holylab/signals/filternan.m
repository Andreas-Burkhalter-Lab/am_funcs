function y = filternan(b,a,x)
% FILTERNAN: filter with an IIR along one dimension, resetting at NaNs
%
% Syntax:
%   y = filternan(b,a,x);
% where
%   b and a are the two sets of IIR filter coefficients (see FILTER);
%   x is the input data (may be multidimensional, must be of type single);
% and
%   y contains the filtered data.
%
% This function attempts to reduce edge transients by calculating the
% filter delays at steady-state; these delays are applied at the initial
% edge, and re-applied any time a "NaN" is encountered in the data.
%
% See also: FILTER.
  
% Copyright 2006 by Timothy E. Holy
  
b = b(:);
a = a(:);
nfilt = max(length(b),length(a));
if (length(b) < nfilt)
  b(nfilt) = 0;
end
if (length(a) < nfilt)
  a(nfilt) = 0;
end
% The following comes from "filtfilt"; see citation therein.
zi = ( eye(nfilt-1) - [-a(2:nfilt) [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
     ( b(2:nfilt) - a(2:nfilt)*b(1) );
% Now do the real work
y = filternan_mex(b,a,x,zi);
