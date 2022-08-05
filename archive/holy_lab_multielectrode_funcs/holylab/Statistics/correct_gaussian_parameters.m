function [N,covarfac,pvalue] = correct_gaussian_parameters(d,n,r2max_over_r2mean)
% correct_gaussian_parameters: estimate the complete Gaussian from partial inflation
% Given a Gaussian centered at zero, that within a radius rmax there were n
% samples drawn, and that the mean square displacement of these points is
% r2mean. The task is to estimate the actual width of the Gaussian and the
% total number of samples drawn.
%
% Syntax:
%   [N,covarfac,pvalue] = correct_gaussian_parameters(d,n,r2max_over_r2mean)
% where
%   d is the dimensionality
%   n is the number of samples drawn up to radius rmax
%   r2max_over_r2mean is the ratio rmax^2/r2mean (r2mean = mean square
%     displacement of the points)
% and
%   N is the estimated total number of points drawn
%   covarfac is the estimated factor by which to multiply the covariance to
%     get the true covariance (i.e., sigma^2 = covarfac*r2mean).
%   pvalue tests the hypothesis that the corrected parameters are
%     indistinguishable from the "original" (technically, with N=n+1). Low
%     pvalue implies that the N,covarfac estimate is more reliable than the
%     original.

% Copyright 2011 by Timothy E. Holy

  func = @(x) -loglikelihood(x,d,n,r2max_over_r2mean);
%   x0 = [1;0.9*d];
%   ops = optimset(@fmincon);
%   ops.Algorithm = 'interior-point';
%   ops.Display = 'none';
%   [x,fval] = fmincon(func,x0,[],[],[],[],[1;0],[inf;d],[],ops);
%   m = x(1); betap = x(2);
  x0 = [0;0];
  ops = optimset(@fminunc);
  ops.LargeScale = 'off';
  ops.Display = 'none';
  [x,fval] = fminunc(func,x0,ops);
  [m,betap] = x2p(x,d);
  N = m+n;
  covarfac = d/betap;
%   if (nargout > 2)
%     % Find the optimal solution with m = 1 (would prefer m=0, but Stirling
%     % isn't so good there...)
%     funcbeta = @(y) func([1;y]);
%     [y1,fval1] = fmincon(funcbeta,0.9*d,[],[],[],[],0,d,[],ops);
%     z = fval1-fval;
%     pvalue = gammainc(z/2,d/2,'upper');
%   end
  
% % Parametrization that is good for fmincon  
% function L = loglikelihood(x,d,n,r2max_over_r2mean)
%   % Joint N,beta form
%   m = x(1);
%   betap = x(2); % betap = beta*R^2
%   logNchoosen = log((1/n+1/m)/2/pi)/2 + (1+n/m+m/n)/(12*(m+n)) + ...
%     n*log(1+m/n) + m*log(1+n/m);
%   L = logNchoosen  + m*log(gammainc(betap*r2max_over_r2mean/2,d/2,'upper')) + ...
%     d*(n+1)*log(betap)/2 - n*betap/2;

% Parametrization that is good for fminunc
function L = loglikelihood(x,d,n,r2max_over_r2mean)
  % Joint N,beta form
  [m,betap] = x2p(x,d);
  logNchoosen = log((1/n+1/m)/2/pi)/2 + (1+n/m+m/n)/(12*(m+n)) + ...
    n*log(1+m/n) + m*log(1+n/m);
  L = logNchoosen  + m*log(gammainc(betap*r2max_over_r2mean/2,d/2,'upper')) + ...
    d*(n+1)*log(betap)/2 - n*betap/2;

function [m,betap] = x2p(x,d)
  m = exp(x(1))+1;
  betap = (atan(x(2))+pi)/(2*pi)*d; % betap = beta*R^2
