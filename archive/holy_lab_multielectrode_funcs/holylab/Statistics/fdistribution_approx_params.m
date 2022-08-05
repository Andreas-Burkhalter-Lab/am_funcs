function [s,n1,n2] = fdistribution_approx_params(mu1,mu2,arg3,type)
% fdistribution_approx_params: calculate parameters for approximate F distribution
% Statistics that involve T^2 variables (or sums of T^2 variables) can be
% exactly or approximately considered to be F-distributed,
%     s*T^2 ~ F(n1,n2).
% This function calculates the parameters s, n1, and n2 from the moments of
% T^2.  See T. E. Holy, "Statistics of the mean using constrained sample
% covariance matrices", 2011 (?).
%
%   [s,n1,n2] = fdistribution_approx_params(mu1,mu2,n1,'2-moment')
% uses the 2-moment approximation, given the mean (mu1) and variance (mu2)
% of T^2. In this case you specify n1 directly (the redundant syntax is to
% preserve compatibility with the 3-moment mode).
%
%   [s,n1,n2] = fdistribution_approx_params(mu1,mu2,mu3,'3-moment')
% uses the 3-moment approximation, additionally given the third central
% moment mu3 (i.e., E((T^2-mu1)^3) where E is expected value).
%
% Note that all three parameters need to be positive and finite, but that
% this is not guaranteed for all inputs, so you need to check the output
% yourself.
%
% See also: fstat, T2_moments.

% Copyright 2011 by Timothy E. Holy

  r1 = mu1^2/mu2;
  switch(type)
    case '3-moment'
      mu3 = arg3;
      r2 = mu1*mu2/mu3;
      n1 = 4*(r1*r2+r1-r2)/(4*r2-r1+1);
      n2 = (4*r1*r2+6*r1-8*r2)/(r1-2*r2);
    case '2-moment'
      n1 = arg3;
      n2 = (2*n1*r1+4*n1-4*r1)/(n1-2*r1);
    otherwise
      error(['Type ' type ' not recognized'])
  end
  s = n2/mu1/(n2-2);
    