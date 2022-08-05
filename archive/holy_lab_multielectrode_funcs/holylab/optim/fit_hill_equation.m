function [results,options] = fit_hill_equation(C,r,sigma,options)
% FIT_HILL_EQUATION: fit measurements to enzyme kinetic model (e.g., Michaelis-Menten)
%
% This function generates fits using a sequence of models of increasing
% complexity:
%   flat: the rate is independent of concentration. This is appropriate if
%     the system is unresponsive or saturated.
%   linear: the rate is linear in concentration. This is appropriate if the
%     system is not near saturation.
%   full: a full fit to the Michaelis-Menten or Hill equation 
%
% The full model used is a first-order process,
%    r = rmax*c/(c + K)
%
% Optionally, you can also fit a Hill exponent n,
%    r = rmax*c^n/(c^n + K^n)
%
% You can also allow an offset,
%    r = r0 + rmax*...
%
% Syntax:
%   results = fit_hill_equation(c,r,sigma)
%   results = fit_hill_equation(C,r,sigma)
%   results = fit_hill_equation(...,options)
%   [results,options] = fit_hill_equation(...)
% where
%   r and sigma are vectors of length n_measurements vectors, containing
%     the mean and s.e.m. of the measured rates for each sample 
%   c is a vector containing the concentration of each sample.
%     Alternatively, to handle multiple compounds, C is a matrix of size
%     n_samples-by-n_components, where C(i,j) is the concentration of the
%     jth component in the ith sample.  In this latter case the full model
%     for the jth sample is
%        r = rmax * ceff^n / (1 + ceff^n)
%     where
%        ceff = sum(cv./Kv)
%     with cv and Kv being the vectors containing the concentrations and
%     binding constants for each component.
%   options is a structure which may have the following fields:
%    fit_n (default false): if true, allows n to be different from 1
%    allow_offset (default false): if true, the baseline (r0) can be
%      nonzero, and is optimized.
%
% On output,
%   results is a structure, with fields "flat", "linear", and "full."  Each
%     of these holds a structure. All 3 have fields "chisq" and "dof",
%     giving the chisq and # of degrees of freedom for the fit to the
%     corresponding model. (You can use these to help select among models.)
%     Also included are parameter values, and their either the inverse
%     covariance matrix ("icov") or, for the "flat" model, the standard
%     error.
%     For the flat model, the predicted firing rate is always rbar.
%     For the linear model, the predicted firing rate is
%         results.linear.r0 + C * results.linear.m(:)
%       given a vector of concentrations c of all the components.
%     For the full model, use hill_equation.
%
% See also: hill_equation, RATE1STORDER, CONCCOMPARE.

% Copyright 2011 by Timothy E. Holy

  %% Argument parsing
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'fit_n',false,'allow_offset',false);
  
  if ~isvector(r)
    error('r must be a vector');
  end
  n_samples = length(r);
  if numel(sigma) ~= n_samples
    error('r and sigma must match');
  end
  r = r(:);
  sigma = sigma(:);
  if any(C(:) < 0)
    error('Concentrations must be non-negative');
  end
  if isvector(C)
    C = C(:);
  end
  if (size(C,1) ~= n_samples)
    error('c must be a vector, or a n_samples-by-n_components matrix');
  end
  % Check for any useless components, with all C zero
  n_components = size(C,2);
  if any(all(C == 0,1)) && n_components>1
    error('Cannot have components whose concentrations are all zero');
  end

  %% Defaults on the bounds on the fit
  cmax = max(C,[],1);
  Cmin = C; Cmin(Cmin == 0) = Inf;
  cmin = min(Cmin,[],1);  % minimum of non-zero values of C
  % upper bounds
  if ~isfield(options,'smax')
    options.smax = struct;
  end
  options.smax = default(options.smax,'rmax',Inf,'K',1000*cmax,...
    'n',10,'r0',Inf);
  % lower bounds
  if ~isfield(options,'smin')
    options.smin = struct;
  end
  options.smin = default(options.smin,'rmax',0,'K',cmin/1000,...
    'n',0.1,'r0',-Inf);
  
  
  %% Initialize the output
  % Order in each model type: parameter values, uncertainty in
  % parameters (icov = inverse covariance), and then fitting quality (chisq/dof)
  results.flat = struct('rbar',NaN,...
    'rbar_err',NaN,...
    'chisq',NaN,'dof',n_samples-1);
  n_var = n_components + options.allow_offset;
  results.linear = struct('m',NaN(1,n_components),'r0',NaN,...
    'parameters','',...
    'icov',NaN(n_var,n_var),...
    'chisq',NaN,'dfof',n_samples-n_var);
  n_var = 1+n_components+options.fit_n+options.allow_offset;
  results.full = struct('rmax',NaN,'K',nan(1,n_components),'n',NaN,'r0',NaN,...
    'parameters','',...
    'icov',NaN(n_var,n_var),...
    'chisq',NaN,'dof',n_samples-n_var);
  if ~options.fit_n
    results.full.n = 1;
  end
  if ~options.allow_offset
    results.linear.r0 = 0;
    results.full.r0 = 0;
  end
  
  if any(sigma == 0)
    % At least one point has zero error bars, so numerical fit is
    % meaningless (and probably the data are weird anyway)
    return
  end
  
  %% Fit to a flat model (appropriate for nonresponsive/saturated)
  % Calculate the sigma-weighted mean response
  Si = 1./sigma.^2;
  results.flat.rbar = sum(r.*Si)/sum(Si);
  % Calculate the confidence limits on rbar
  results.flat.rbar_err = 1/sqrt(sum(Si));
  % See whether the variation across samples is significant
  results.flat.chisq = sum((r-results.flat.rbar).^2.*Si);
  
  if all(C == 0)
    return
  end
  
  %% Fit to a linear model (appropriate if not saturated)
  % m = rmax/K must be nonnegative, but r0 can have any sign
  rsc = r./sigma;
  Csc = bsxfun(@rdivide,C,sigma);
  lb = zeros(1,n_components);
  ub = inf(1,n_components);
  if options.allow_offset
    % Allow a baseline different from zero. We handle this with one
    % additional row to the matrix
    Csc = [Csc, 1./sigma];
    lb(end+1) = -inf;
    ub(end+1) = inf;
  end
  results.linear.m = lsqlin(Csc,rsc,[],[],[],[],lb,ub,[],optimset('Display','off'))';
  results.linear.icov = Csc'*Csc;
  dr = rsc - Csc*results.linear.m(:);
  results.linear.chisq = sum(dr.^2);
  results.linear.parameters = 'm';
  if options.allow_offset
    % The final coordinate is the baseline response, peel it off
    results.linear.r0 = results.linear.m(end);
    results.linear.m = results.linear.m(1:end-1);
    results.linear.parameters = {'m','r0'};
  end
  
  %% Fit to the full Hill equation
  % Create initial guesses:
  % rmax
  s0.rmax = max(r);
  % K
  % This is tricky because don't want K == inf, so in cases where 1/K is
  % tiny, replace K by some factor times cmax
  Ki = results.linear.m / s0.rmax;  % inverse K
  Kthresh = 10*cmax;
  s0.K = min(1./Ki,Kthresh);
  % n
  s0.n = 1;
  % r0
  s0.r0 = results.linear.r0;
  % Choose which parameters are to be optimized
  fn = {'rmax','K'};
  if options.fit_n
    fn = [fn {'n'}];
  end
  if options.allow_offset
    fn = [fn {'r0'}];
  end
  results.full.parameters = fn;
  % Set up the functions converting the structure to a parameter vector and
  % vice versa
  extract_func = @(s) extract_fields(s,fn{:},struct('reshape',false));
  [p0,fields,field_shape,sbase] = extract_func(s0);
  fill_func = @(p) fill_fields(fields,field_shape,p,sbase);
  lb = extract_func(options.smin);
  ub = extract_func(options.smax);
  % Define the optimization function
  myfun = @(s) hill_penalty(s,C,r,sigma);
  opt_func = @(p) optimize_struct_wrapper(p,myfun,extract_func, ...
    fill_func);
  % Perform the nonlinear least squares optimization
  [p,results.full.chisq,~,~,~,~,J] = lsqnonlin(opt_func,p0,lb,ub,optimset('Jacobian','on','Display','off','DerivativeCheck','off'));
  results.full.icov = J'*J;
  s = fill_func(p);
  results.full = copyfields(s,fn,results.full);
  

function [resid,Js] = hill_penalty(s,C,r,sigma)
  ceff = sum(bsxfun(@rdivide,C,s.K),2);
  ceffn = ceff.^s.n;
  ceffn1 = ceffn ./ (1+ceffn);
  rpred = s.r0 + s.rmax * ceffn1;
  resid = (rpred - r)./sigma;
  if (nargout > 1)
    Js.rmax = ceffn1./sigma;
    Js.r0 = 1./sigma;
    denom = (1+1./ceffn).^2 .* ceffn .* sigma;
    Js.n = s.rmax*log(ceff)./denom;
    Js.n(ceff == 0) = 0;
    Js.K = C .* ((s.n*s.rmax./denom./ceff) * (-1./s.K.^2));
    Js.K(ceff == 0,:) = 0;
  end

