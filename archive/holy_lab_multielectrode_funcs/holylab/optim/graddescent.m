function [p,err] = graddescent(func,p0,options)
% GRADDESCENT: a gradient descent minimization algorithm
%
% Syntax:
%   p = graddescent(func,p0,options)
% where
%   func is the function handle to your function to be minimized (which
%     has a syntax [val,grad] = func(p), where val is the value and grad
%     is the gradient with respect to the parameters);
%   p0 is your initial guess for the parameters
%   options is a structure which may have the following fields:
%     iter_max (default 100): the maximum number of parameter tries;
%     mu (default 1): the initial factor multiplying the gradient in
%       determining the step size;
%     mu_decfactor (default 5): the factor by which to decrease mu on
%       failed steps;
%     mu_incfactor (default 2): the factor by which to increase mu on
%       successful steps;
%     tol (default 1e-6): the fractional/absolute change in the function
%       value needed for convergence;
%     step_min: if present, uses the step size as the alternative
%       criterion for convergence. Each parameter must change by less
%       than its corresponding entry in the vector step_min.
%     goodstep_fcn: a function handle to execute whenever a good step is
%       found. The syntax is f(p), where p is the current parameter state.
%
% See also: FMINUNC.
  
% Copyright 2006 by Timothy E. Holy
	
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'iter_max')
    options.iter_max = 100;
  end
  if ~isfield(options,'mu')
    options.mu = 1;
  end
  if ~isfield(options,'mu_decfactor')
    options.mu_decfactor = 5;
  end
  if ~isfield(options,'mu_incfactor')
    options.mu_incfactor = 2;
  end
  if ~isfield(options,'tol')
    options.tol = 1e-6;
  end
  options = default(options,'linesearch_plot',false,'Display',false);
  
  if options.linesearch_plot
    figure
    hax = gca;
  end
  p = p0;
  mu = options.mu;
  iter = 0;
  isdone = false;
  err = [];  % history of error
  [v,g] = func(p);
  while ~isdone
    if (options.linesearch_plot & ishandle(hax))
      alpha = linspace(-5*mu,5*mu,11);
      val = nan(size(alpha));
      for alphaIter = 1:length(alpha)
        ptest = p - alpha(alphaIter)*g;
        val(alphaIter) = func(ptest);
      end
      axes(hax)
      plot(alpha,val)
      shg
      pause
    end
    % Try going along the gradient until the function value diminishes.
    % Note that we aren't really trying a line minimization here; it
    % seems that line minimization might slow down progress compared with
    % adaptive step size approaches.
    % The issue comes down to the ratio in expense between computing
    % function values and gradients: if the gradient is not much worse
    % than the function value, it seems one should compute more
    % gradients and fewer function values, since you get more information
    % from the gradient.
    ptest = p - mu*g;
    g2 = g(:)'*g(:);
    vNew = func(ptest);
    if (vNew >= v)
      % Bad step, diminish mu
      mu = mu/options.mu_decfactor;
      % Should we quit?
      if isfield(options,'step_min')
        if all(mu*abs(g) < options.step_min)
          isdone = true;
        end
      else
        if (mu*g2 < (abs(v)+1)*options.tol)
          isdone = true;
        end
      end
    else
      % Good step
      p = ptest;
      err(end+1) = v;
      [v,g] = func(p);
      if options.Display
        fprintf('iter %d: err %g, mu = %g\n',iter,vNew,mu);
      end
      if isfield(options,'goodstep_fcn')
        options.goodstep_fcn(p);
      end
      mu = mu*options.mu_incfactor;
    end
    iter = iter+1;
    if (iter > options.iter_max)
      isdone = true;
    end
  end
  if (iter > options.iter_max)
    warning('graddescent:premature','convergence not achieved');
  end
  err(end+1) = v;
  %iter
