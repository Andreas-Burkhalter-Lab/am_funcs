function [x,fvalret] = mgopt(func,data,x0,options)
% MGOPT: multigrid optimization of nonlinear functions
%
% Multigrid is a technique for handling problems with large numbers of
% variables. It is best known for the solution of partial differential
% equations, but it is also applicable to optimization problems in which
% the variables live on a "grid." The central idea is to use coarse-grid
% solutions, which can be computed more quickly, to improve fine-scale
% solutions.
%
% This implements the MG/OPT algorithm described in the following paper:
%   "A multigrid approach to discretized optimization problems",
%   Steven G. Nash, Optimization Methods and Software, 2000.
% It does not constrain the line search, and for some problems this turns
% out to be a problem; the gradient correction ends up dominating the
% problem and it yields poor results.
%
% Syntax:
%   [x,fval] = mgopt(func,data,x0,options)
% where
%   func: a function handle with syntax
%         [val,grad] = func(x,data1,data2,...)
%     is the function whose value you want to minimize with respect to x.
%     The gradient (grad) must be in the shape of the grid (not a vector as
%     with MATLAB's routines). "func" can also depend upon data whose size
%     also varies with the grid (data1, data2, etc).  Any grid-independent
%     variables can be specified using anonymous functions.
%   data: an n_levels-by-n_data cell array of arrays, where the n_levels
%     are the data restricted to increasingly coarse grids (i.e., data{1,1}
%     is the first data array at the finest resolution, and data{5,1} would
%     be the same array restricted four times to coarser grids).
%     MG_PREPARE_DATA is a convenience function that can help you prepare
%     coarser versions of your data. Note that it's not required that you
%     generate coarser data by restriction; if there's an analytic model
%     that would be better for coarse-grid data, you can use it.
%   x0: supply an initial guess or leave empty (in which case it will start
%     solving on the coarsest grid with an all-zeros guess).  You can
%     supply a guess of a size equal to any of the grids (finest or any
%     coarse-grid approximation).
%   options: a structure array with the following fields:
%     N0, N1: the number of pre- and post-polishing steps at each phase
%       (default 1)
%     tol (default 1e-6): the fractional change in return value required
%       for convergence
%     n_cycles: if set, simply executes the given number of multigrid
%       V-cycles (tol is ignored in this case)
% and
%   x is the optimized solution
%   fval is a vector giving the function value on each V-cycle.
%
% Example: smoothing a vector c by minimizing
%    E[u] = int dx (u - c)^2 + lambda * int dx (grad u)^2 
% ("int" means integral). This is equivalent to solving a Helmhotz
% equation. Higher lambda means greater smoothing, since it sets the
% penalty on the gradient of u.
%
% You could generate a vector c in the following way:
%    npts = 1025; t = linspace(0,6*pi,npts); c = sin(t) + randn(1,npts);
%    lambda = 1000;
% and the coarse-grid versions this way:
%   data = {c,lambda};
%   while (length(data{end,1}) > 24)
%     data{end+1,1} = array_restrict(data{end,1},[false true]);
%     data{end,2} = data{end-1,2}/4;  % the correct scaling for this problem
%   end
%
% and set up "func" in this way:
%    [val,grad] = helmholtz_func(u,c,lambda)
%       du = diff(u);
%       val = sum((u-c).^2) + lambda * sum(du.^2);
%       grad = -2*lambda*lap1d(u) + 2*(u - c);
% (lap1d is presumed to be a 1-dimensional laplacian function).
%
% Then you can solve the problem with
%   [u,fval] = mgopt(@(u,c,lambda) helmholtz_func(u,c,lambda),data);
% You should see that it exits after the 2nd V-cycle.
%
% For comparison, try
%   [ucg,fvalcg] = conjgrad(@(u) helmholtz_func(u,c,lambda),zeros(size(c)),struct('iter_max',1000));
% It will exit in comparable time, but with slightly higher function value
% even after 130 iterations. In general, multigrid and conjugate-gradient
% have their own strengths, and you may want to experiment to find the best
% choice for your particular problem. Multigrid's advantages tend to
% increase for multidimensional and/or large-scale problems.
%
% See also: MG_PREPARE_DATA, ARRAY_RESTRICT, ARRAY_PROLONG.

% Copyright 2009 by Timothy E. Holy
  
  %% Parse the arguments
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'N0',1,'N1',1,'tol',1e-6,'display_cg',false,'display_mg',true,'protect_v',true,'mu_start',1,'f0',0);
  if ~iscell(data)
    error('Data must be supplied as an n_levels-by-n cell array');
  end
  [n_levels,n_data] = size(data);
  if (n_levels == 1)
    if (n_data > 1)
      warning('multigrid:dataFormat','Assuming you meant to supply data as a column cell array');
      data = data(:);
      n_levels = n_data;
      n_data = 1;
    else
      error('Only one level supplied.');
    end
  end
  have_starting_guess = nargin > 2 && ~isempty(x0);
  if have_starting_guess
    n_dims = ndims(x0);
  elseif isfield(options,'n_dims')
    n_dims = options.n_dims;
  else
    n_dims = zeros(1,n_data);
    for dataIndex = 1:n_data
      n_dims(dataIndex) = ndims(data{1,dataIndex});
    end
    n_dims = min(n_dims);
  end
  mode = 'tolerance';
  if isfield(options,'n_cycles')
    mode = 'n_cycles';
  end

  %% Calculate the grid size & restriction schedule at each level
  grid_size = zeros(n_levels,n_dims);
  schedule = false(n_levels-1,n_dims);
  for lvlIndex = 1:n_levels
    sz = size(data{lvlIndex,1});
    sz(end+1:n_dims) = 1;
    grid_size(lvlIndex,:) = sz(1:n_dims);
    if (lvlIndex > 1)
      schedule(lvlIndex-1,:) = grid_size(lvlIndex,:) ~= ...
        grid_size(lvlIndex-1,:);
    end
  end
  
  %% Assign the starting level
  if have_starting_guess
    level = [];
    for lvlIndex = 1:n_levels
      if isequal(size(x0),grid_size(lvlIndex,:))
        level = lvlIndex;
        break
      end
    end
    if isempty(level)
      error('Starting guess does not correspond to any grid size');
    end
  else
    x0 = zeros(grid_size(end,:));
    level = n_levels;
  end
  
  %% Do the minimization, testing for convergence
  x = x0;
  if strcmp(mode,'tolerance')
    [x,fval] = mg_vcycle(x,level);
    fvalOld = fval;
    fvalret = fval;
    isdone = false;
    while ~isdone
      [x,fval] = mg_vcycle(x,1);
      isdone = fvalOld - fval < options.tol * (abs(fvalOld) + abs(fval));
      fvalOld = fval;
      fvalret(end+1) = fval; %#ok<AGROW>
    end
  else
    fvalret = zeros(1,options.n_cycles);
    for cycleIndex = 1:options.n_cycles
      [x,fvalret(cycleIndex)] = mg_vcycle(x,level);
      level = 1;
    end
  end
  
  %% Utility functions
  % The function value, corrected for the gradient
  function [fval,grad] = mg_func_corrected(x,dat,v,xbase,omega)
    [fval,grad] = func(x,dat{:});
    dx = x - xbase;
    grad = grad - v + dx .* omega;
    fval = fval - sum(x(:) .* v(:)) + sum(dx(:).^2 .* omega(:))/2;
  end

  % Linesearch on a corrected function
  function [mu,fval] = mg_linesearch(x,dx,mu,dat,v,xbase,omega,fval)
    deltax = x - xbase;
    linfunc = @(mu) func(x + mu * dx,dat{:}) ...
	      - sum((x(:) + mu * dx(:)) .* v(:)) ...
        + sum((deltax(:) + mu * dx(:)).^2 .* omega(:))/2;
    [mu,fval] = linmin(linfunc,struct('startval',fval,'mu',mu));
  end

  % "Improve" by doing one step of gradient-descent
  function [x,fval,mu] = mg_improve(x,dat,v,xbase,omega,mu)
    [fval,grad] = mg_func_corrected(x,dat,v,xbase,omega);
    [mu,fval] = mg_linesearch(x,-grad,mu,dat,v,xbase,omega,fval);
    x = x - mu*grad;
  end

  %% The core multigrid functions
  function [x,fval] = mg_vcycle(x,level)
    z = zeros(size(x));
    [x,fval] = mg_recurse(x,level,z,z);
    while (level > 1)
      % When starting with an initial guess on a coarse grid, we need to
      % work our way up to finer grids
      level = level - 1;
      x = array_prolong(x,grid_size(level,:));
      z = zeros(size(x));
      [x,fval] = mg_recurse(x,level,z,z);
    end
  end

  function [x,fval] = mg_recurse(x,level,v,omega)
    xbase = x;
    if (level == n_levels)
      % We're on the lowest level, solve the problem by
      % conjugate-gradient
      cgfunc = @(x) mg_func_corrected(x,data(n_levels,:),v,xbase,omega);
      [x,fval] = conjgrad(cgfunc,x,struct('iter_max',1000,'Display',options.display_cg,'mu',options.mu_start));
      if options.display_mg
        fprintf('%d CG iterations\n',length(fval)-1);
      end
    else
      % Do a coarse-grid correction
      % Step 1: polish the current guess at the current level
      mu = options.mu_start;
      for polishIndex = 1:options.N0
        [x,fval,mu] = mg_improve(x,data(level,:),v,xbase,omega,mu);
      end
      % Step 2: calculate the term needed for correcting the gradient at the
      % coarser level
      [fval,gradh] = mg_func_corrected(x,data(level,:),v,xbase,omega);  % need to get the gradient
      flag = schedule(level,:);
      xH = array_restrict(x,flag);
      z = zeros(size(xH));
      [valH,gradH] = mg_func_corrected(xH,data(level+1,:),z,z,z);
      gradhr = array_restrict(gradh,flag);
      vH = gradH - gradhr;
      %vH = 0*vH;
      if options.protect_v
        vH(gradH == 0) = 0;  % in case some variables irrelevant on coarse scale
      end
      omegaOld = array_restrict(omega,flag);
      omegaNew = vH.^2/(abs(fval)+abs(valH)+options.f0);
      omegaH = omegaOld+omegaNew;
      % Step 3: solve on the coarser grid (recursive)
      xHnew = mg_recurse(xH,level+1,vH,omegaH);
      % Step 4: bring the correction up to the finer grid
      dx = array_prolong(xHnew-xH,grid_size(level,:));
      % Step 5: optimize the use of the coarse-grid correction
      [mu,fval] = mg_linesearch(x,dx,1,data(level,:),v,xbase,omega,fval);
      if options.display_mg
        fprintf('level %d: mu = %g\n',level,mu);
      end
      x = x + mu*dx;
      % Step 6: perform fine-grid polishing
      mu = options.mu_start;
      for polishIndex = 1:options.N1
        [x,fval,mu] = mg_improve(x,data(level,:),v,xbase,omega,mu);
      end
    end
  end

end


% Note: might want to bound the change in x by +-gamma*e, where e is the
% all-ones vector and gamma a scalar,
%   gamma = max{|vH|,|restrict(gradh)|,|gradH|}
% MODEL PROBLEMS FOR THE MULTIGRID OPTIMIZATION OF
%SYSTEMS GOVERNED BY DIFFERENTIAL EQUATIONS
%  ROBERT MICHAEL LEWIS AND STEPHEN G. NASH  2005
% Can easily implement this at the level of line search

