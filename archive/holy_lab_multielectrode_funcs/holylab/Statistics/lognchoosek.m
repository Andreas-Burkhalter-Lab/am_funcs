function ret = lognchoosek(n,k)
% lognchoosek: Stirling's approximation for nchoosek
% choose(n,k) can get big quickly. This function uses Stirling's
% approximation to give an approximate value of log(choose(n,k)).
%
% Syntax:
%   ret = lognchoosek(n,k)
%
% See also: nchoosek.

% Copyright 2011 by Timothy E. Holy

  m = n-k;
  ret = log((1/n+1/m)/2/pi)/2 + (1+n/m+m/n)/(12*(m+n)) + ...
    n*log(1+m/n) + m*log(1+n/m);
