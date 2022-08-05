function [p,err] = conjgrad(func,p0,options)
% CONJGRAD: a conjugate-gradient descent minimization algorithm
%
% Syntax:
%   p = conjgrad(func,p0,options)
% where
%   func is the function handle to your function to be minimized. If you
%     are doing unconstrained minimization, this has a syntax
%         [val,grad] = func(p),
%     where val is the value and grad is the gradient with respect to the
%     parameters. If you are doing constrained minimization (using a
%     barrier penalty), then this should have a syntax
%         [val,grad] = func(p,barrier_mu)
%   p0 is your initial guess for the parameters
%   options is a structure which may have the following fields:
%     iter_max (default 100): the maximum number of parameter tries;
%     mu (default 1): the initial factor multiplying the gradient in
%       determining the step size;
%     tol (default 1e-6): the fractional/absolute change in the function
%       value needed for convergence;
%     preconditioner: if supplied, preconditioned nonlinear conjugate
%       gradient is used. The preconditioner should be a symmetric
%       positivie definite approximation of the Hessian that reduces the
%       condition number of the true Hessian. The syntax is 
%           s = options.preconditioner(p,g)
%       where p is the current parameter value, g is the current gradient,
%       and s is the altered gradient. Default is s = g (which you
%       get by not supplying the preconditioner).
%       If you're using a barrier function (see below), the syntax is
%           s = options.preconditioner(p,barrier_mu,g)
%     geodesic_fcn: provides the opportunity to search with equality
%       constraints. The syntax is
%           pstar = options.geodesic_fcn(p,h,mu),
%       where p is the current position, h is the search direction (which
%       may or may not be tangent to the manifold at p), and mu is the
%       amplitude. Under unconstrained conditions, the candidate position
%       would be
%            p + mu*h,
%       but this function allows you to modify this location to some
%       altered position pstar which is on the manifold. The only
%       requirement is that the mapping from mu to pstar needs to be
%       (locally) one-to-one and twice continuously differentiable.
%     barrier (default false): if true, performing constrained minimization
%       using a barrier function (inequality constraints). In this case,
%       the next two parameters are also used:
%         barrier_mu: the initial value of the barrier penalty mu.
%         barrier_mu_step (default 0.9): the new value of barrier_mu on
%           each step, relative to the previous value:
%             barrier_mu_new = barrier_mu_step * barrier_mu_old
%     goodstep_fcn: a function handle to execute whenever a good step is
%       found. The syntax is f(p), where p is the current parameter state.
%       To save parameter values at each step of the iteration, see
%       PHISTORY.
%     Display: if true, prints convergence information.
%     DisplayParams: if true, prints parameter values on each accepted step
%
% [p,fval] = conjgrad(...)
%   also returns the function value on each iteration of the algorithm.
%
% Based on the exposition in J. Shewchuk, "An introduction to the conjugate
% gradient method without the agonizing pain."
%
% See also: FMINUNC, GRADDESCENT, PHISTORY.
  
% Copyright 2008 by Timothy E. Holy
	
  if (nargin < 3)
    options = struct;
  end
  options = default(options,...
    'iter_max',100,...
    'mu',1,...
    'tol',1e-6,...
    'preconditioner',[],...
    'geodesic_fcn',[],...
    'barrier',false,...
    'barrier_mu_step',0.9,...
    'barrier_contrib_small',true,...
    'linesearch_plot',false,...
    'Display',false,...
    'DisplayParams',false,...
    'DisplayParamsIndex',1:length(p0));
  
  use_preconditioner = ~isempty(options.preconditioner);
  use_geodesic = ~isempty(options.geodesic_fcn);
  if options.linesearch_plot
    figure
    hax = gca;
  end
  p = p0;
  mu = options.mu;
  iter = 0;
  isdone = false;
  err = [];  % history of error
  
  % Compute the initial function value and search direction
  if ~options.barrier
    [v,g] = func(p);
  else
    barrier_mu = options.barrier_mu;
    [v,g] = func(p,barrier_mu);
  end
  if (all(g(:) == 0)) % can not check inequity of g(:) with 0.0 in this way
    p = p0;
    err = v;
    return
  end
  if isinf(v) || isnan(v)
    error('Must have finite starting value of function');
  end
  if options.Display
    fprintf('iter    0: value = %8.6e \n',v);
  end
  if options.DisplayParams
    fprintf('%g ',p(options.DisplayParamsIndex));
    fprintf('\n');
  end
  if use_preconditioner
    if ~options.barrier
      s = options.preconditioner(p,g);
    else
      s = options.preconditioner(p,mu,g);
    end
  else
    s = g;
  end
  h = s;  % the search direction

  % Begin iterations
  deltaNew = g(:)'*s(:);  % this is g^2 if not using preconditioning
  while ~isdone
    if (options.linesearch_plot && ishandle(hax))
      % Display the data used for line search (useful for debugging)
      alpha = linspace(-5*mu,5*mu,11);
      val = nan(size(alpha));
      for alphaIter = 1:length(alpha)
        ptest = p - alpha(alphaIter)*h;
        val(alphaIter) = func(ptest);
      end
      plot(hax,alpha,val)
      shg
      pause
    end
    % Do the line minimization along -h
    % This differs from Shewchuk, in that it should provide guaranteed
    % convergence (we explicitly require a minimum of the function value
    % rather than a zero of the gradient along the line).
    if use_geodesic
      if options.barrier
        fmu = @(mu) func(options.geodesic_fcn(p,-h,mu),barrier_mu);
      else
        fmu = @(mu) func(options.geodesic_fcn(p,-h,mu));
      end
    else
      if ~options.barrier
        fmu = @(mu) func(p - mu*h);
      else
        fmu = @(mu) func(p - mu*h,barrier_mu);
      end
    end
    [mutriple,fmutriple] = linbrackgd(fmu,v,-sum(g(:).*h(:)),mu);
    [mu,vfinal] = brentmin(fmu,mutriple,fmutriple,struct('TolFun',max(options.tol/10,sqrt(eps))));
    if use_geodesic
      p = options.geodesic_fcn(p,-h,mu);
    else
      p = p - mu*h;
    end
    % Compute the value and gradient at the new point
    if ~options.barrier % if isnan(vNew) or isfinite(vNew) true, then start next
      % calculation, instead of continue this loop
      [vNew,gNew] = func(p);
    else
      barrier_mu = barrier_mu * options.barrier_mu_step;
      [vNew,gNew] = func(p,barrier_mu);
    end
    if options.Display
      if ~options.barrier
        fprintf('iter %4d: value = %8.6e, mu = %8.6e \n',iter+1,vNew,mu);
      else
        fprintf('iter %4d: value = %8.6e, mu = %8.6e, barrier_mu = %8.6e\n', ...
          iter+1,vNew,mu,barrier_mu);
      end
    end
    if options.DisplayParams
      fprintf('%g ',p(options.DisplayParamsIndex));
      fprintf('\n');
    end
    if isfield(options,'goodstep_fcn')
      options.goodstep_fcn(p);
    end
    % Compute the new search direction (conjugate to previous directions)
    deltaOld = deltaNew;
    if (deltaOld == 0) % can not compare float number in this way
      isdone = true;
    end
    deltaMid = gNew(:)'*s(:);
    if use_preconditioner
      if ~options.barrier
        s = options.preconditioner(p,gNew);
      else
        s = options.preconditioner(p,mu,gNew);
      end
    else
      s = gNew;
    end
    deltaNew = gNew(:)'*s(:);
    beta = (deltaNew-deltaMid)/deltaOld;
    if (beta < 0 || mod(iter+1,numel(p)) == 0)
      h = gNew;
    else
      h = gNew + beta*h;
      if (sum(gNew(:).*h(:)) < 0)
        h = gNew;
      end
    end
    % Test for convergence
    % This criterion also differs from Shewchuck, it uses only the function
    % value (the square gradient may not even make sense from a units
    % perspective)
    if (v - vfinal < options.tol*(abs(v)+abs(vfinal))) % use abs(v - vfinal)?
      if options.barrier && options.barrier_contrib_small
        % If using a barrier function, make sure the barrier contributes
        % very little
        v0 = func(p,0);
        if (abs(vfinal-v0) < options.tol*(abs(v0)+abs(vfinal)))
          isdone = true;
        end
      else
        isdone = true;
      end
    end
    if (mu == 0)
      isdone = true;
    end
    % Prepare for the next iteration
    err(end+1) = v; %#ok<AGROW>, store convergence history
    v = vNew;
    g = gNew;
    iter = iter+1;
    if (iter >= options.iter_max)
      isdone = true;
    end
  end
  if (iter >= options.iter_max)
    warning('conjgrad:premature','convergence not achieved');
  end
  err(end+1) = v;
