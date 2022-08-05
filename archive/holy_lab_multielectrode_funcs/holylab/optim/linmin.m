function [mu,fval] = linmin(f,options)
% LINMIN: one-dimensional minimization along a line
% If you have a multidimensional function f, a starting point x0, and a
% downhill direction v, then the function
%    fmu = @(mu) f(x0 + mu*v)
% is a one-dimensional function of mu that decreases (at least for a
% while) for mu > 0. LINMIN finds the minimum of fmu with respect to mu.
%
% Syntax:
%   mu = linmin(fmu)
%   [mu,fval] = linmin(fmu,options)
% where
%   fmu is the anonymous function of mu
%   options is a structure which may have the following fields:
%     mu (default 1): first value of mu to try
%     mu_max: if supplied, controls the largest value of mu permitted in a
%       step
%     startval: if supplied, it uses this value for fmu(0) (thus skipping
%       the calculation, useful if calculation is slow and you already
%       know the value at mu=0)
%     bracketonly (default false): if true, returns the central mu in a
%       bracketing operation, but does not optimize mu
%     itermax (default 10): if the initial mu is too large, mu will be
%       decreased 10-fold to try to find a value smaller than fmu(0), up to
%       itermax times
% and
%   mu is the value of mu at a minimum of fmu
%   fval is fmu(mu).
%
% See also: LINDIP, LINDEC, GRADDESCENT.
  
% Copyright 2007 by Timothy E. Holy
  
  % Calls to optimset can take a lot of time, so do this only once
  persistent ops
  if isempty(ops)
    ops = optimset(@fminbnd);
  end
  
  if (nargin < 2)
    options = struct;
  end
  if (~isfield(options,'mu') || isempty(options.mu))
    options.mu = 1;
  end
  if isfield(options,'startval') && ~isempty(options.startval)
    startval = options.startval;
  else
    startval = f(0);
  end
  if (isinf(startval) || isnan(startval))
    error('linmin cannot work with Inf or NaN values');
  end
  options = default(options,'bracketonly',false,'itermax',10,'mu_max',inf);
  if isfield(options,'TolFun')
    ops.TolX = options.TolX;
  end
  
%   if isfield(options,'mu_max')
%     ops.TolX = 1e-4;
%     [mu,fval] = fminbnd(f,0,options.mu_max,ops);
%     fval
%     return
%   end
  
  lastmu = 0;
  lastlastmu = 0;
  lastval = startval;
  lastlastval = startval;
  mu = options.mu;
  
  % Our first task is to bracket the minimum
  % Check to make sure we haven't exceeded the bounds of the function
  tempval = f(mu);
  while (isinf(tempval) || isnan(tempval))
    mu = mu/10;
    tempval = f(mu);
  end

  % Progress along positive mu until we start to rise again
  while (tempval < lastval && mu < options.mu_max)
    lastlastmu = lastmu;
    lastmu = mu;
    lastlastval = lastval;
    lastval = tempval;
    mu = mu*10;
    tempval = f(mu);
  end
  if (mu > options.mu_max)
    mu = options.mu_max;
  end
  
  % If the best mu so far is 0, we need to backtrack until we find a place
  % that improves the guess---otherwise, we might get stuck making no
  % improvement at all
  if (lastmu == 0)
    lastmu = mu/10;
    lastval = f(lastmu);
    iter = 0;
    while (lastval > startval && iter < options.itermax)
      mu = lastmu;
      lastmu = mu/10;
      lastval = f(lastmu);
      iter = iter+1;
    end
    if (lastval > startval)
      % No progress was made
      mu = 0;
      fval = startval;
      return
    end
  end
  
  % Now we have the minimum bracketed between lastlastmu and mu
  if options.bracketonly
    mu = lastmu;
    fval = lastval;
  else
    % Find a local minimum in mu
    ops.TolX = 1e-3*abs(mu-lastlastmu);
    [mu,fval] = fminbnd(f,lastlastmu,mu,ops);
    % Since fminbnd never evaluates the endpoints, we could get into some
    % trouble here---check to see if the solution is worse than the
    % endpoints
    while (fval > lastlastval)
      mu = lastlastmu + 0.5*(mu-lastlastmu);
      [mu,fval] = fminbnd(f,lastlastmu,mu,ops);
    end      
  end
  if (fval > startval)
    error('Something is wrong!'); % don't error, it might be roundoff?
  end
    