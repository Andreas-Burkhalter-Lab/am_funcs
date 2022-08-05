function [mu,fval] = lindec(f,options)
% LINDEC: function value decrement in one dimension (loose minimization)
% If you have a multidimensional function f, a starting point x0, and a
% downhill direction v, then the function
%    fmu = @(mu) f(x0 + mu*v)
% is a one-dimensional function of mu that decreases (at least for a
% while) for mu > 0. LINDEC simply finds a value of mu so that
%       fmu(mu) < f(0).
% This is an alternative to finding the true minimum with respect to mu
% that may lead to faster code in some circumstances.
%
% The algorithm is very simple: start with some guess for mu. If
% fmu(mu) < f(0), you're done; if not, keep making mu smaller until the
% condition is satisfied.
%
% Code that calls this might want to increase the starting guess for mu
% (say by a factor of 2) for the next iteration, so that over time a good
% setting for mu is found adaptively.
%
% Syntax:
%   mu = lindec(fmu)
%   [mu,fval] = lindec(fmu,options)
% where
%   fmu is the anonymous function of mu
%   options is a structure which may have the following fields:
%     mu (default 1): first value of mu to try
%     startval: if supplied, it uses this value for fmu(0) (thus skipping
%       the calculation, useful if calculation is slow and you already
%       know the value at mu=0)
%     max_iter (default 40): the maximum number of times this will try to
%       decrease mu before giving up and throwing an error
%     mu_decrease_factor (default 5): factor by which to decrease mu when
%       fmu(mu) > f(0)
% and
%   mu is the value of mu at a minimum of fmu
%   fval is fmu(mu).
%
% See also: LINDIP, LINMIN, GRADDESCENT.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  mu_set = false;
  if (~isfield(options,'mu') || isempty(options.mu))
    mu_set = true;
    options.mu = 1;
  end
  options = default(options,'max_iter',40);
  options = default(options,'mu_decrease_factor',5);
  if (isfield(options,'startval') && ~isempty(options.startval))
    startval = options.startval;
  else
    startval = f(0);
  end
  mu = options.mu;
  
  fval = f(mu);
  if mu_set
    % This is our first call. Keep boosting mu until the function value
    % increases
    while (fval < startval)
      mu = 10*mu;
      fval = f(mu);
    end
  end
  iter = 0;
  while (fval > startval && iter < options.max_iter)
    mu = mu / options.mu_decrease_factor;
    fval = f(mu);
    iter = iter+1;
  end
  if (fval > startval)
    error('optim:maxiter',...
	  ['Maximum number of iterations exceeded without decreasing' ...
	   ' function value']);
  end

  