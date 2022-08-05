function [n,mu,R,y] = fit_truncated_gaussian(x,domainfunc,params,n,mu,R)
  % fit_truncated_gaussian: fit a Gaussian density to points over an incomplete domain
  %
  % Suppose we have a set of points sampled from a distribution, but that
  % these points are also required to fall within a particular region to be
  % observed.  This function allows you to fit the observed data points to
  % a Gaussian defined over the allowed region.
  % The mathematics behind this is described in the supplementary
  % information to T. E. Holy, "some paper," 2012.
  %
  % Syntax:
  %   [n,mu,R,y] = fit_truncated_gaussian(x,domainfunc)
  %   [n,mu,R,y] = fit_truncated_gaussian(x,domainfunc,n0,mu0,R0)
  % where
  %   x is a d-by-N matrix of points that fall within the specified domain
  %   domainfunc is a function, domainfunc(y) is true if the d-by-1 vector
  %     falls within the domain, and false otherwise
  %   params has fields:
  %     covarianceModel: 'isotropic', 'diagonal', or 'full'
  %     avoidFull (default false): set true if it's impractical to invert a
  %       d-by-d matrix
  %     nmax (default inf): maximum value for n
  %   n0, mu0, R0 (optional) are an initial guess for the number of points,
  %     the center, and the inverse square bandwidth of the Gaussian.
  % On output, the Gaussian is specified as
  %     f(x) = |detR|/(2*pi)^(d/2) exp(-sum((R*(x-mu)).^2)/2)
  % The output y is a set of points sampled from this Gaussian that are all
  % within the specified domain.

  % Copyright 2012 by Timothy E. Holy
  
  %% Initialize
  params = default(params,'avoidFull',false,'nmax',inf);
  if (params.avoidFull && strcmp(params.covarianceModel,'full'))
    error('Cannot set avoidFull=true and also use full covariance model');
  end
  [d,nx] = size(x);
  thiseps = sqrt(eps(class(x)));  % numerical precision threshold
  if nargin < 4
    n = nx;
  end
  if nargin < 5
    mu = mean(x,2);
  end
  dx = bsxfun(@minus,x,mu);
  if nargin < 6
    Cx = (dx*dx')/n; % not n-1 because ML yields n
    R = chol(Cx);
    R = inv(R)';
  end
  
%   %% Iteratively improve the estimate---Expectation-Maximization method
%   % Commented out because it doesn't work very well! Very slow
%   % convergence, but more importantly there must be something wrong because
%   % the log-likelihood seems to decrease rather than increase.
%   Lold = inf;
%   L = inf;
%   while (~isfinite(L) || L < Lold)
%     Lold = L;
%     % Generate sufficient random points that we get nx points within the
%     % domain
%     sigma = 1/sqrt(beta);
%     y = bsxfun(@plus,sigma*randn(d,ceil(n)),mu);
%     keep = domainfunc(y);
%     ck = cumsum(keep);
%     n_new = find(ck >= nx,1,'first');
%     while isempty(n_new)
%       n_extra = (length(ck)/ck(end) - 1)*nx;
%       y_extra = bsxfun(@plus,sigma*randn(d,ceil(n_extra+5*sqrt(n_extra))),mu);
%       keep_extra = domainfunc(y_extra);
%       ck_extra = cumsum(keep_extra)+ck(end);
%       keep = [keep keep_extra];
%       ck = [ck ck_extra];
%       y = [y y_extra];
%       n_new = find(ck >= nx,1,'first');
%     end
%     n = n_new;
%     F = nx/n;
%     keep = keep(1:n);
%     y = y(:,1:n);
%     % Replace the points within the domain with x
%     y(:,keep) = x;
%     % Compute the optimum parameters
%     mu = mean(y,2);
%     dy = bsxfun(@minus,y,mu);
%     sdy2 = sum(dy(:).^2);
%     beta = d*n/sdy2;
%     fprintf('mu:');
%     fprintf(' %g',mu);
%     fprintf(', beta: %g\n',beta);
%     % Compute the log-likelihood
%     L = n*d*log(beta/2/pi)/2 - beta*sdy2/2;
%   end
  

  %% Iteratively improve the estimate---Newton's method
  % The logic to maximize the log-likelihood is a bit complicated, because
  % we want to ensure global convergence, the Monte Carlo introduces
  % statistical uncertainty, and we want to avoid code-duplication. The
  % main issue is that on a Newton step, the log-likelihood should
  % increase; if it doesn't that means either (1) we took too large of a
  % step, or (2) we're so close to the peak that the value is dominated by
  % Monte Carlo sampling error. In the first case we want to take a smaller
  % step, and on the second we want to terminate. We tell when we've
  % transitioned from 1 to 2 by the heuristic that once we accept the full
  % Newton step, we're likely in the domain of convergence (this is not
  % fool-proof and may need to be re-evaluated).
  Lold = -inf; % the likelihood using our best-yet parameters
  converging = false;  % this will be true once we accept the full Newton step
  fullnewton = true;
  first = true;        % is this the initial pass? (If so, we should not turn on "converging")
  hfig = figure;
  while true
    L = nan;     % the current likelihood (nan if not yet determined)
    iter_inner = 0;
    while isnan(L) || (L < Lold)
      % Perform the Monte Carlo
      % Draw the same number of points from the current Gaussian estimate
      iter_inner = iter_inner+1;
      sigma = inv(R);
      [y,n] = randn_indomain(nx,n,mu,sigma,domainfunc,params.nmax,0.1);
%       plot(x(1,:),x(2,:),'b.');
%       line(y(1,:),y(2,:),'LineStyle','none','Marker','.','Color','r')
%       n
%       mu
%       sigma
%       pause
      if (size(y,2) == nx)
        % Compute the mean log-likelihood (ignoring constants like 2*pi)
        F = nx/n;
        dx = bsxfun(@minus,x,mu);
        Rdx = R*dx;
        Lx = log(abs(det(R))) - mean(sum(Rdx.^2,1))/2;
        L = Lx - log(F);
      end
      if (L < Lold || size(y,2) < nx)
        if converging && size(y,2) == nx
          % We're within the basin of convergence, this likely indicates
          % sampling uncertainty. There is little point in continuing.
          close(hfig)
          return
        else
          if first
            % We haven't even successfully calculated L yet, but the
            % initial draw did not yield a sufficient "hit rate". Make
            % sigma smaller
            R = 2*R;
          else
            % We're not in the basin of convergence, so we likely took too
            % large a step. Cut the step size by a factor of 2.
            fullnewton = false;  % next time it will not be full Newton step
            deltanu = deltanu/2;
            deltaR = deltaR/2;
            R = Rold + deltaR;
            deltamu = R\deltanu;
            mu = muold + deltamu;
          end
        end
      else
        % It was the full Newton step, so signal that from now on we can
        % assume we're within the basin of convergence
        if fullnewton && ~first
          converging = true;
        end
      end
      if (iter_inner > 20)
        return  % after this many tries, couldn't do better, so we must be p<0.05 of the optimum
      end
    end
    %% Draw a fresh data set from the same model
    % By doing this before we store L, it ensures that we don't select
    % "increasingly-lucky" data sets with successive iterations.
    [y,n] = randn_indomain(nx,n,mu,sigma,domainfunc,inf);
    % Recompute the likelihood
    F = nx/n;
    L = Lx - log(F);
    % Prepare for the next iteration by storing old parameters and setting
    % state
    muold = mu;
    Rold = R;
    Lold = L;
    first = false;  % we are no longer in the initial calculation of the likelihood
    fullnewton = true;
    % Compute the gradient and a semi-definite approximation to the Hessian
    dy = bsxfun(@minus,y,mu);
    Rdy = R*dy;
    S = bsxfun(@times,reshape(Rdy,[d 1 nx]),reshape(dy,[1 d nx]));
    S = reshape(S,[d^2 nx]);
    meanRdy = mean(Rdy,2);
    meanS = mean(S,2);
    gnu = mean(Rdx,2) - meanRdy;
    gr = -(Rdx*dx' - Rdy*dy')/nx;
    negHnunu = (Rdy*Rdy')/nx - meanRdy*meanRdy';
    negHnur = -(Rdy*S')/nx + meanRdy*meanS';
    negHrr = (S*S')/nx - meanS*meanS';
    % Solve for the update
    g = [gnu; gr(:)];
    H = [negHnunu, negHnur; negHnur' negHrr];
    % Because R can be multiplied by an orthogonal matrix without consequence,
    % H is guaranteed to be poorly conditioned. So the inversion needs to
    % be done carefully.
    [U,S] = svd(H);
    s = diag(S);
    keep = s > thiseps;  % this should discard d*(d-1)/2 components
    sinv = zeros(size(s));
    sinv(keep) = 1./s(keep);
    deltaparam = U*diag(sinv)*U'*g;
    deltanu = deltaparam(1:d);
    deltaR = reshape(deltaparam(d+1:end),[d d]);
    % Prevent the diagonals of R from going negative---this isn't
    % mathematically inadmissable, but it does represent a sharp change in
    % where the simulated points end up and the fraction that stay within
    % the permitted region.
    Rdiag = diag(Rold);
    dRdiag = diag(deltaR);
    negflag = Rdiag+dRdiag < 0;
    if any(negflag)
      % Permit a step to half of the current value
      alpha = -Rdiag(negflag)./dRdiag(negflag)/2;
      alpha = min(alpha);
      deltaR = alpha*deltaR;
      deltanu = alpha*deltanu;
      fullnewton = false;
    end
    % Update the parameters
    R = Rold + deltaR;
    deltamu = R\deltanu;
    mu = muold + deltamu;
  end
end

%% Draw target number of in-domain points from Gaussian
% Draw a specified # of points from the Gaussian distribution that lie
% within the domain. The target # of points is nx. mu and sigma are the
% parameters of the Gaussian, the input n is the current guess for the # of
% points we need to draw. On output, y contains the in-domain points and n
% is the total number of draws that were required.
% Thresh is optional; if supplied, this terminates early if we've kept
% fewer than thresh*nx points using our initial guess n.
function [y,n] = randn_indomain(nx,n,mu,sigma,domainfunc,nmax,thresh)
  d = length(mu);
  n = ceil(n);
  y = bsxfun(@plus,sigma*randn(d,n),mu);
  keep = domainfunc(y);
  y = y(:,keep);
  % Did we select a sufficient number of the first round?
  if (size(y,2) >= nx)
    y = y(:,1:nx);
    ck = cumsum(keep);
    n = find(ck >= nx,1,'first');
    return
  end
  % Check to see if we're generating fewer in-domain points than
  % minimally-acceptable, in case we have set the parameters poorly.
  if (nargin > 6)
    if (size(y,2) < nx*thresh)
      return  % terminate early
    end
  end
  % Add more points until we get enough, keeping track of the total number
  % sampled
  while (size(y,2) < nx)
    n_extra = (n/size(y,2) - 1)*nx; % expected # of extra needed
    n_extra = n_extra+5*sqrt(n_extra);  % an overabundance so we likely can finish quickly
    n_extra = min(n_extra,nmax-n);
    if (n_extra == 0)
      break
    end
    y_extra = bsxfun(@plus,sigma*randn(d,ceil(n_extra)),mu);
    keep = domainfunc(y_extra);
    ck = cumsum(keep);
    n_criterion = find(size(y,2) + ck >= nx,1,'first');
    if isempty(n_criterion)
      n_criterion = n_extra;
    else
      keep(n_criterion+1:end) = false;
    end
    n = n+n_criterion;
    y = [y y_extra(:,keep)];   %#ok<AGROW>
  end
end
