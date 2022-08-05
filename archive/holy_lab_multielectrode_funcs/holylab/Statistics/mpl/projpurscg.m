function w = projpurscg(x,n,dp,options)
% PROJPURSCG: projection pursuit by conjugate gradient
% By default this does PPEK (Projection Pursuit with Exponential Kernels)
% for reasons of speed.
% 
% Syntax:
%   w = projpurscg(x,n,dp,options)
% where
%   x contains the data points, and is of size d-by-Nunique (d is the dimensionality);
%   n contains the multiplicity of each point, and is of size 1-by-Nunique;
%   dp is the dimension of the projection (e.g., if you want to reduce to 2
%     dimensions, set dp = 2);
% and
%   w is a dp-by-d matrix containing the projection directions, one per
%     row.  Note that these directions will not be unit vectors (they're
%     scaled by the optimal magnitude for 1-d kernel-density estimation),
%     nor will they necessarily be orthogonal. They will, however, be
%     conjugate with respect to the covariance matrix of the data
%     points (i.e., w*Cx*w' is diagonal to within roundoff).
% You can control aspects of this calculation through the structure
% "options," which may have the following fields:
%     sequential (default true): if true, finds the dimensions one at a
%       time, so you're doing PPEK. Warning: if you set this to be false
%       and dp>1, you're using an O(N^2) rather than an O(N) algorithm! (It
%       would be possible to implement a simultaneous 1d projection
%       algorithm, but this code doesn't do that.) Be prepared for it to be
%       slow. But the advantage is the projection pursuit is being carried
%       out for higher-dimensional projections, which in theory might allow
%       you to discover patterns that you would otherwise miss.
%     start (default 'random'): defines the starting value for each
%       projection. 'random' does what it suggests; 'pca' chooses each
%       initial projection direction from a principal components analysis
%       of the data; you can also supply a numeric input, in which case the
%       initial w = options.start.
%     wprev: if present, it will assume that wprev represents
%       previously-chosen projection directions. The found directions will
%       be conjugate (Cx-orthogonal) to wprev (so that w*Cx*wprev' will
%       be 0). This allows you to find the projection directions
%       incrementally.
%     oversmooth (default 0.5): controls the amount of oversmoothing
%       applied in seeking new directions. oversmooth of 0 means no
%       oversmoothing (the width of the kernels is chosen optimally); in
%       general the width of the kernels is (1+oversmooth)*(optimal width);
%
%     tol (default 1e-4): the tolerance on the angle change in two
%       successive iterations before convergence is declared;
%     minsep (default 1e-4): when doing cross-validation to optimize the
%       length of w, "too-close" points are consolidated. "Too close" is
%       judged as points which are less than 1e-4*(total range of projected
%       data). In other words, the data are placed on a grid with 1/minsep
%       points in it.
%     itermax (default 100+4*d): the maximum number of iterations
%       acceptable;
%     mag_maxiter (default 20): the maximum number of iterations for which
%       the magnitude will be adjusted. Once the magnitude is fixed, the
%       direction will still be optimized until the convergence criterion
%       is achieved.
%
%     verbose: off (no text output), iter (a small amount at each
%       iteration), or full (quite a lot at each iteration)
%     Cx: supply the correlation matrix of the data x if you want to skip
%       the time required to calculate it (or if d is big, and you don't
%       have the storage for Cx: supply pre-whitened data and Cx as a
%       sparse diagonal matrix).
%
% See also: MPL_OPTW.

% Copyright 2006 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end

  if ~isfield(options,'sequential')
    options.sequential = 1;
  end
  if ~isfield(options,'start')
    options.start = 'random';
  end
  if ~isfield(options,'verbose')
    options.verbose = 'iter';
  end
  
  [d,Np] = size(x);
  N = sum(n);

  % Compute the covariance matrix of the data. This will be useful for
  % computing the gradient of the gaussian-fit-to-projection, as well
  % as for providing an initial guess for w by PCA.
  % Note we do it in this way to make use of the multiplicity
  % information, rather than just calling svd on x.
  % Alternatively, it might be smarter to simply insist that the user
  % whiten the inputs? In high dimensions, Cx is going to be big...
  if ~isfield(options,'Cx')
    nd = spdiags(n',0,Np,Np);  % multiplicity info
    xm = mean(nd*x')';
    dx = x - repmat(xm,1,Np);
    Cx = (dx*nd*dx')/(N-1); % Note that if this is too high-dimensional, pre-whiten & supply as sparse
  else
    Cx = options.Cx;
  end
  % Check for previously-set directions
  if isfield(options,'wprev')
    wprev = options.wprev;
  else
    wprev = zeros(0,d);
  end
  nprev = size(wprev,1);
  % Pick starting directions
  if isnumeric(options.start)
    w0 = options.start;
  elseif ischar(options.start)
    switch options.start
      case 'random'
        w0 = randn(dp,d);
      case 'pca'
        [U,S,V] = svd(Cx);  % do PCA
        dS = diag(S);
        w0 = diag(1./sqrt(dS((1:dp)+nprev)))*V(1:d,(1:dp)+nprev)'; % w is dp-by-d
      otherwise
        error(['Starting option ' options.start ' not recognized']);
    end
  else
    options.start
    error(['Starting option not recognized']);
  end

  if (options.sequential || dp == 1)
    % Sequential 1-d projections
    w = [wprev;zeros(size(w0))];
    for i = 1:dp
      w0cur = w0(i,:);
      % Minimize the action subject to the constraint of being orthogonal to
      % the previously-found directions
      [wnew,negen] = ppwork(w0cur,w,x,n,Cx,options);
      w(i+nprev,:) = wnew;
      negentropy(i) = negen;  % Don't return this, it's not right because of oversmoothing, etc.
    end
  else
    % Multidimensions, simultaneous
    [w,negentropy] = ppwork(w0,zeros(size(w0)),x,n,Cx,options);
  end
  
function [w,negen] = ppwork(w0,w,x,n,Cx,options)
  if (nargin < 6)
    options = struct;
  end
  if ~isfield(options,'oversmooth')
    options.oversmooth = 0.5;
  end
  N = sum(n);
  [dp,d] = size(w0);
  
  % Compute Cx-unit vectors parallel to the previous proj. directions
  %wnorm = sqrt(sum(w.^2,2));
  wnorm = sqrt(diag(w*Cx*w'));  % by construction from previous iterations, w*Cx*w' is diagonal
  isz = (wnorm == 0);
  wnorm(isz) = 1;  % avoid divide-by-zero
  what = w .* repmat(1./wnorm,1,d);  % These are Cx-unit vectors (or 0)
  
  % Start so that projection is orthogonal to prev. directions
  w0 = w0 - (w0*Cx*what')*what;
  
  % Optimize the overall length of w0 at the starting point
  X = w0*x;
  mplops = struct('kernel_only',1,'mode','diag');
  if ~strcmp(options.verbose,'full')
    mplops.display = 'none';
  end
  [Q,wscale] = mpl_optw(X,n,mplops);
  w0 = wscale*w0;
  
  if ~isfield(options,'tol')
    tol = 1e-4;
  else
    tol = options.tol;
  end
  if ~isfield(options,'minsep')
    options.minsep = 1e-4;
  end
  if ~isfield(options,'itermax')
    itermax = 100+4*d;
  else
    itermax = options.itermax;
  end
  if ~isfield(options,'mag_maxiter')
    options.mag_maxiter = 20;
  end
  iter = 0;
  reset_counter = 0;
  diranglestep = 0.1*pi;  % a tenth of the way around half-circle
  magstep = 1/N;
  isdone = 0;
  smoothfac = 1+options.oversmooth;
  
  [val2,grad2] = pp_mnegentropy(w0/smoothfac,x,n,Cx,what);
  r = -grad2;  % residual
  % precondition here?
  conj = r;
  deltanew = trace(r*conj');
  while ~isdone
    val0 = val2;
    grad0 = grad2;
    % Approximately minimize along the direction conj from w0
    dirops = struct('normalize',1,'val0',val0,'grad0',grad0,'verbose',options.verbose);
    diranglestepOld = diranglestep;
    [diranglestep,w0tmp,val1] = linmincg(@(w) pp_mnegentropy(w,x,n,Cx,what),...
                                  w0/smoothfac,conj,2*diranglestep,dirops);
    w0 = w0tmp*smoothfac;                                
    % Before we complete the conjugate-gradient on the direction, we also
    % need to optimize the overall length of the projection vector(s).
    % Note: this step is the most annoying numerically, due to the issues
    % with cross-validation. Consider turning this off after a certain
    % number of iterations?
    if (iter < options.mag_maxiter)
      % Compute the projected data
      xp = w0*x;
      % Aggregate too-close points
      [xp,np] = uniquemult(xp',n,options.minsep);
      xp = xp';
% $$$       if (dp == 1)
% $$$         [xp,sort_order] = sort(xp);
% $$$         np = n(sort_order);
% $$$         % Consolidate near-duplicates by rounding to an integer lattice
% $$$         rng = diff(xp([1 end]));
% $$$         xp = round((xp-xp(1))/rng/options.minsep)+1;
% $$$         [clabel,nlabel] = agglabel(xp);
% $$$         xp = find(nlabel);
% $$$         np = nlabel(xp);
% $$$         xp = xp*rng*options.minsep;
% $$$         magops.sorted = 1;
% $$$       else
% $$$         [xp,sort_order] = sortrows(xp');
% $$$         np = n(sort_order);
% $$$         % Consolidate near-duplicates here?
% $$$         magops.sorted = 0;
% $$$       end
      mag0 = ones(1,dp);
      [val1p,grad1p] = pp_mloglikeCV(mag0,xp,np);
      magops = struct('normalize',0,...
        'val0',val1p,'grad0',grad1p,'verbose',options.verbose);
      [magstep,mag] = linmincg(@(w) pp_mloglikeCV(w,xp,np),...
        mag0,-grad1p,2*magstep,magops);
      % Multiply the projection direction(s) by these magnitudes
      w0 = w0.*repmat(mag',1,size(w0,2));
    else
      mag = 1;
      magstep = 0;
    end
    if strcmp(options.verbose,'full')
      fprintf('iter %d, val0 %g, val1 %g, diranglestep %g, magstep %g, mag ',iter,val0,val1,diranglestep,magstep);
      fprintf('%g ',mag);
      fprintf(', w0:\n');
    elseif strcmp(options.verbose,'iter')
      fprintf('.');
    end
    % OK, now we can get back to the conjugate-gradient analysis
    [val2,grad2] = pp_mnegentropy(w0/smoothfac,x,n,Cx,what);
    r = -grad2;  % residual
    % precondition here?
    deltaold = deltanew;
    deltamid = trace(r*conj');
    deltanew = trace(r*r');
    beta = (deltanew-deltamid)/deltaold;
    reset_counter = reset_counter+1;
    if (reset_counter == d || beta <= 0)
      conj = r;
      reset_counter = 0;
      % beta <= 0 seem associated with pathologically small steps
      if (beta <= 0)
        diranglestep = max([diranglestep diranglestepOld]);
      end
      if strcmp(options.verbose,'full')
        fprintf('reset\n');
      end
    else
      conj = conj - ((conj*w0')/(w0*w0'))*w0; % Subtract component || to w0
      % (since w0 has changed by looping around a circle, what was once
      % perpendicular to w0 will no longer be)
      conj = r + beta*conj;
    end
    if (r*conj' < 0 || diranglestep < tol)
      % Approximate line minimizations were not close enough to true
      % minimum, so we have to re-start
      conj = r;
      reset_counter = 0;
      if strcmp(options.verbose,'full')
        fprintf('reset2\n');
      end
    end
    iter = iter+1;
    %isdone = (abs(val2-val0) + abs(val1-val0) < (abs(val0)+tol)*tol || ...
    %          iter > itermax);
    % Require that the diranglestep is less than tolerance for two
    % successive iterations
    isdone = ((diranglestep < tol && ...
      diranglestepOld < tol && abs(mag-1) < tol)  || iter > itermax);
  end
  if ~strcmp(options.verbose,'none')
    fprintf(' (%d iterations)\n',iter);
  end
  if (iter > itermax)
    warning('projpurscg failed to converge')
  end
  w = w0;
  negen = -val2;

function [val,grad] = pp_mnegentropy(w,x,n,Cx,what)
% Compute -negentropy (plusentropy?), and its gradient
% We want to minimize -negentropy
  [dp,d] = size(w);
  CY = w*Cx*w';
  N = sum(n);
  if (dp == 1)
    % One dimensional projection
    if (nargout > 1)
      [S,gradS] = llcv1(w,x,n);
      gradS = gradS';
      gradSg = (N/CY)*(w*Cx);
    else
      S = llcv1(w,x,n);
    end
    Sg = (N/2)*log(CY) + N*(1 + log(2*pi))/2;
  else
    % Multidimensional projection
    Y = w*x;
    sd = sqrdist(Y,Y);
    G = exp(-sd/2)/(2*pi)^(d/2);
    N = sum(n);
    Q = (G*n')/N;
    S = -n*log(abs(Q)); % log likelihood, w/o det(w)
    % Now do the gaussian with the covariance of the projected data
    CY = w*Cx*w';
    Sg = (N/2)*log(det(CY)) + N*dp*(1 + log(2*pi))/2;
    if (nargout > 1)
      gradS = gradwllcv(x,Y,n,G,Q,'full'); % gradient of the above
      gradSg = N*(CY\(w*Cx));
    end
  end
  val = S-Sg;
  if (nargout > 1)
    grad = gradS - gradSg;
    if (nargin > 4)
      % Subtract off component Cx-parallel to previous directions
      grad = grad - (grad*Cx*what')*what;
      % Subtract off component that is || to w (which changes the magnitude)
      grad = grad - ((grad*w')/(w*w'))*w;
    end
  end

function [val,grad] = pp_mloglikeCV(wdiag,x,n)
% Compute -loglikelihood, and its gradient
  dp = length(wdiag);
  wdiag = abs(wdiag);
  N = sum(n);
  if (dp == 1)
    % One dimensional projection
    if (nargout > 1)
      [val,grad] = llcv1(wdiag,x,n,struct('CV',1));
      grad = grad - N/wdiag;
    else
      val = llcv1(wdiag,x,n,struct('CV',1));
    end
    val = val - N*log(wdiag);
  else
    % Multidimensional projection
    X = repmat(wdiag',1,size(x,2)).*x;
    sd = sqrdist(X,X);
    G = exp(-sd/2)/(2*pi)^(dp/2);
    Gcv = G;
    for i = 1:size(G,1)
      Gcv(i,i) = 0;
    end
    Pcv = (Gcv*n')/N;
    val = -n*log(abs(Pcv)) - N*sum(log(wdiag));
    if (nargout > 1)
      grad = gradwllcv(x,X,n,G,Pcv,'diag'); % gradient of the above
      grad = grad' - N./wdiag;
    end
  end
  
  
function [step_opt,wmin,funcval] = linmincg(func,w0,dw,stepinc,options)
% Approximately minimize along the direction dw from w0
% (once the minimum is bracketed, we'll just interpolate its position
% and call it good enough)
  val0 = options.val0;
  grad0 = options.grad0;
  if isfield(options,'iterlm_max')
    iterlm_max = options.iterlm_max;
  else
    iterlm_max = 100;
  end
  if (options.normalize)
    % We're going along a circle around the origin. Thus, steps can be
    % expressed as an angle, and in fact we know the maximum angle worth
    % trying (pi). To make the parametrization easy, make dw have the same
    % length as w0. Then our test points will be
    %     cos(theta)*w0 + sin(theta)*dw
    stepinc_max = pi;
    dw = dw*sqrt((w0*w0')/(dw*dw'));
  else
    stepinc_max = inf; % We're doing line search. We can go as far as we want.
    dwtest = dw; % The tangent along the line will always be dw
  end
  slopedec = trace(dw*grad0');  % Slope along line while decreasing
  if ~(slopedec < 0)
    warning('slopedec is not negative');
    step_opt = 0;
    wmin = w0;
    funcval = val0;
  end
  stepdec = 0;
  iterlm = 0;
  isdone = 0;
  while ~isdone
    % Test a new point, and see if we've "rounded the bend" and started
    % going uphill again
    if options.normalize
      % We're going around the circle
      wtest = cos(stepinc)*w0 + sin(stepinc)*dw;
      dwtest = -sin(stepinc)*w0 + cos(stepinc)*dw;
    else
      % We're doing a real line search
      wtest = w0 + stepinc*dw;
    end
    [val1,grad1] = func(wtest);
    slopeinc = trace(dwtest*grad1');  % Find place with increasing slope
    iterlm = iterlm+1;
    % Note we also need to check the function value: we could have
    % overshot the minimum & a maximum, and be back to decreasing slope
    % even if the function value is larger than where we started
    isdone = ~(slopeinc < 0 && val1 < val0 && ...
      iterlm < iterlm_max && stepinc < stepinc_max);
    if (~isdone && val1 < val0 && slopeinc < 0)
      % OK, we went further downhill. Keep track of this, because it's
      % useful progress towards the minimum.
      stepdec = stepinc;
      slopedec = slopeinc;
      stepinc = 2*stepinc;
      stepinc = min(stepinc,stepinc_max);
    end
    % Now get back to trying to find where it starts going uphill
  end
  % We've either exhausted our patience (iterlm == iterlm_max), or we've
  % found the turning point. Either way, we want to use the information
  % we've got to find a point that's lower in value than when we
  % started---and if we have found the turning point, we can use the slope
  % information to interpolate the zero-crossing, i.e., the true
  % minimum. As long as that interpolated position reduces the function
  % value, we'll settle for it, rather than zooming in on the exact minimum
  % in this direction (to reduce the number of function evaluations).
  % Here's the linearly-interpolated position of the minimum, as a
  % proportion of the distance between stepdec & stepinc:
  frac = -slopedec/(abs(slopeinc)-slopedec);
  % Let's see if it's really lower (maybe we _way_ overshot, esp if our
  % initial guess for stepinc was lousy)
  isdone = 0;
  funceval = iterlm+1;
  iterlm = 0;
  while ~isdone
    steptest = (1-frac)*stepdec + frac*stepinc;
    if options.normalize
      wtest = cos(steptest)*w0 + sin(steptest)*dw;
    else
      wtest = w0 + steptest*dw;
    end
    val1 = func(wtest);
    iterlm = iterlm+1;
    % Did we succeed in finding a better point?
    isdone = (val1 < val0 || iterlm > iterlm_max);
    if ~isdone
      % Nope, try something closer to the last point we found where the
      % function value was decreasing (which might well be the start
      % position)
      frac = frac/2;
    end
  end
  if (iterlm > iterlm_max && val1 >= val0)
    % We never did find a better point. This suggests we might already
    % have been at the minimum. Return the starting position.
    steptest = 0;
    wtest = w0;
    val1 = val0;
    warning('Line search failed to find a lower point')
  end
  if strcmp(options.verbose,'full')
    fprintf('Function evals: %d (%d before check)\n',funceval+iterlm,funceval);
  end
  step_opt = steptest;   % We did it!
  wmin = wtest;
  funcval = val1;
