function K = randbindconst(distparms,sz)
% RANDBINDCONST: generate random binding constants from a tilted-log
% distribution
%  
% K = randbindconst(distparms,sz)
% where
%   distparms: [log10Kmin log10Kmax rho], rho is the tilt
%   sz is the size 2-vector

% Generate the binding constants. Assume the distribution
% of log10K is linear over some region from log10Kmin to log10Kmax,
% and zero outside that region
% For x = (log10K-log10Kmin)/(log10Kmax-log10Kmin),
%   p(x) = 2 * (1/(1+rho) * (1-x) + 1/(1+1/rho) * x).
% (Note rho = p(1)/p(0).)
% To draw from this distribution, let y be a uniform deviate.
% Then
%   x = (sqrt(1+(rho^2-1)*y) - 1)/(rho - 1).
  dlK = distparms(2)-distparms(1);
  rho = distparms(3);
  y = rand(sz);
  if (rho == 1)
    x = y;
  else
    x = (sqrt(1+(rho^2-1)*y) - 1)/(rho - 1);
  end
  K = 10.^(distparms(1) + dlK*x);	% These are the binding constants
