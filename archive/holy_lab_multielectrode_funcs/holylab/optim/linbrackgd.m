function [x,fx] = linbrackgd(f,f0,fp0,mu)
% LINBRACKGD: find bracketing triple in gradient-descent 1-d minimization
%
% Approaches like conjugate gradient perform 1-d minimization on a
% function of the following form:
%    f(mu) = F(p - mu*h),
% where F is the multdimensional objective function, p is the base point,
% and h is the search direction.  Furthermore, one also knows the
% function value f(0) = F(p).  This routine is designed to take advantage
% of this information to increase efficiency.  It assumes the minimum
% occurs for some mu>0 (which it should, given the definition of f).
%
% Syntax:
%   [x,fx] = linbrackgd(f,f0,fp0,mu0)
% where
%   f is the 1-dimensional function to be minimized
%   f0 is the function value at mu = 0
%   fp0 is the function derivative at mu = 0; given the expression for f
%     above, you should supply this as fp0 = -sum(g(:).*h(:)), where g is
%     the gradient of F at p.
%   mu0 is a candidate step
% and
%   x, fx contain bracketing-triple information as described in LINBRACK.
%
% See also: LINBRACK, FMINBND.
  
% Copyright 2010 by Timothy E. Holy
  
  iter_max = 20;
  if (fp0 >= 0)
    error('The derivative at mu=0 must be negative');
  end

  fmu = f(mu);
  iter = 0;
  while ((isnan(fmu) || isinf(fmu)) && mu > 0 && iter < iter_max)
    mu = mu/10;
    fmu = f(mu);
    iter = iter+1;
  end
  if (iter >= iter_max)
    warning('No value of mu could be found that yielded finite function values');
    x = 0;
    fx = f0;
    return
  end
  % Use parabolic interpolation to estimate the location of the minimum,
  % given the supplied data
  a = (fmu - mu*fp0 - f0)/mu^2;  % coefficient of the 2nd-order term
  if (a > 0)
    % The function appears to be convex over the interval [0 mu]
    mustar = -fp0/(2*a);
    if (fmu < f0)
      % We can safely use linbrack, because it will search in the
      % direction mu > 0
      [x,fx] = linbrack(f,[0 mu mustar],[f0 fmu]);
    else
      % We have an upper bound for mu, but we need to find a point in the
      % middle that is lower than either end. Note that if we got here,
      % mustar < mu.
      fmustar = f(mustar);
      while (fmustar > f0 && mu > 0 && iter < iter_max)
        mu = mustar;
        fmu = fmustar;
        a = (fmu - mu*fp0 - f0)/mu^2;
        mustar = -fp0/(2*a);  % Note that this will bisect or better
        fmustar = f(mustar);
        iter = iter+1;
      end
      x = [0 mustar mu];
      fx = [f0 fmustar fmu];
    end
  else
    % The function is not convex.  We can just use interval
    % magnification to search.  Note fmu < f0 since fp0 < 0 and a < 0,
    % so linbrack will not search in the mu<0 direction
    [x,fx] = linbrack(f,[0 mu],[f0 fmu]);
  end
end
