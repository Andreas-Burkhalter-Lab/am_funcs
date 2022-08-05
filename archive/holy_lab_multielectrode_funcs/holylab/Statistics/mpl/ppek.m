function [w,dw2tot] = ppek(x,n,options)
% PPEK: Projection Pursuit with Exponential Kernels
% This function performs symmetric projection pursuit, finding the set of
% projection directions which maximizes the sum of the negentropies along
% the individual 1-d projections.  The density is estimated with an
% exponential kernel, which allows for an efficient O(N) algorithm.
%
% Syntax:
%   w = ppek(x,n)
%   w = ppek(x,n,options)
%   [w,dw2tot] = ppek(...)
% where
%   x is a d-by-M matrix of points in d-dimensions, each column a
%     different data point;
%   n is a 1-by-M vector containing the multiplicities of each point;
% and
%   w is a dp-by-d matrix, where each row is a projection direction (see
%     below for how you specify dp, the number of directions);
%   dw2tot is a vector which contains the total sum-of-squares change
%     in w on each iteration.
% The calculation can be controlled in some important ways with options,
% a structure which may have the following fields:
%   dp: the number of projection directions to calculate, up to d (the
%     dimensionality of the data).  Typically, this defaults to
%     d. However, if you specify 'start' (see below) as a matrix, and
%     don't supply this field, then dp is chosen as the number of rows of
%     the starting w.
%   start (default 'random'): either a vector/matrix containing the
%     initial value of w, or a string which specifies how to calculate
%     the starting guess. 'random' chooses gaussian random numbers, 'pca'
%     choses the principal components of x for the initial value.  If you
%     supply a matrix for start, but also provide options.dp > number of
%     rows of options.start, then the remaining projection directions are
%     initialized using 'random'. This provides a mechanism for asking
%     for one more projection direction---but note that in general, the
%     previous directions will change because this is a symmetric
%     algorithm.
%   maxcomponents (default d): controls whether dimensionality reduction
%     (through PCA) is performed before projecting.  Only the projections
%     with the largest maxcomponents eigenvalues are retained.
%   setmag (default false): if true, optimize the length of each w before
%     starting the iterative process of improving both the direction and
%     the magnitude.
%   oversmooth (default 0.5): the degree to which the density should be
%     oversmoothed (compared to the current estimate of the optimum) for
%     the purposes of optimizing direction. Set to 0 for no
%     oversmoothing; 0.5 corresponds to a smoothing length scale 1.5
%     times larger than optimal.
%   outlier_frac (default 0): helps with handling outliers, which can
%     distort the optimal W. Specifies the fraction of points to be
%     treated as outliers, and thus not contribute to the
%     cross-validated log-likelihood. 0 retains all points; 0.05 would
%     discard the 5% least-likely points.
%   random_magnitude: controls whether downhill steps are always taken,
%     or whether the Newton step is computed, but the step taken is then
%     a gaussian random vector centered on the Newton step.  If
%     random_magnitude is empty, then the algorithm will always move
%     downhill (the Newton step is attempted, but a smaller step will be
%     taken if the Newton step increases the action).  If
%     random_magnitude is specified as a vector, this vector determines
%     the schedule for the standard deviation of the random step around
%     the Newton step.  Default: [] if starting with pca or a numeric w0,
%     but has value 4 for a number of trials (minimum of 100 and 2*d),
%     and then linearly decreases to 0 in the next equivalent number of
%     trials.
%   itermax: the maximum number of iterations before giving up, if
%     convergence has not been achieved.
%   display (default 'iter'): controlls the amount of text output. The
%     choices are, from most output to least output, 'full', 'iter', and
%     'none'.
%
% See also: PROJPURSCG.
  
% Copyright 2006 by Timothy E. Holy

% Develop
% Add oversmoothing on the direction?
% Separate out whitening step
% Allow for adding new direction(s)?

% Notes: with the Hessian being approximated as proportional to the
% identity, this is basically a gradient descent, except we have a good
% idea of how large to make the step size.
% Since the Hessian won't really be diagonal, perhaps one could do even
% better with a conjugate-gradient type algorithm?
% Or maybe this is the wrong path: better to go towards stochastic gradient
% descent and make it less vulnerable to local minima? Or perhaps both
% could be provided in different phases.
  
  [d,Np] = size(x);
  if (nargin < 3)
    options = struct;
  end
  options = ppek_default_options(options,d);

  N = sum(n);
  if (d > N)
    warning('Dimensionality greater than the number of points!')
  end
  
  % Check for a numeric start
  if isnumeric(options.start)
    w0 = options.start;
    [dp,w0d] = size(w0);
    if (w0d ~= d)
      error('The sizes of w0 and x do not match');
    end
    if (dp > d)
      error(['This algorithm can''t give you more projection directions than ' ...
        'coordinates.']);
    end
  end
    
  % Center x on zero and then give it unit covariance
  nd = spdiags(n',0,Np,Np);  % multiplicity info
  xm = mean(nd*x')';
  x = x - repmat(xm,1,Np);
  Cx = (x*nd*x')/(N-1);
  [E,D] = eig(Cx);
  options.dp = min([options.dp options.maxcomponents]);
  indxkeep = d:-1:d-options.maxcomponents+1; % Largest eigenvalue is last, make it first
  E = E(:,indxkeep);
  rootDdiag = sqrt(diag(D));
  D = diag(rootDdiag(indxkeep));
  Mwhite = D\E';
  Mwhiteinv = E*D;
  
  % Whiten the data
  x = (Mwhite*x)';   % Note that we've switched x to being N-by-d
  
  % Whiten the starting value, or generate one if need be
  mdp = min([options.dp options.maxcomponents]);
  if isnumeric(options.start)
    w = D*E'*w0';   % Note now the dirs are the columns
    if (options.dp > dp)
      w = [w, randn(options.maxcomponents,options.dp-dp)];
    end
  elseif ischar(options.start)
    switch options.start
      case 'random'
        w = randn(options.maxcomponents,mdp);
      case 'pca'
        w = [eye(mdp); zeros(options.maxcomponents-mdp,mdp)];
      otherwise
        error(['options.start = ' options.start ' not recognized'])
    end
  end
  
  w = w/sqrtm(w'*w);   % To make them orthonormal
  
  if ~isreal(x)
    error('Whitened data points are not real!');
  end
  if ~isreal(w)
    error('Whitened directions are not real!');
  end
  
  if options.setmag
    mplops = struct('kernel_only',1,'mode','diag');
    if ~strcmp(options.display,'full')
      mplops.display = 'none';
    end
    for i = 1:mdp
      X = x*w(:,i);
      [ptmp,wscale] = mpl_optw(X',n,mplops);
      w(:,i) = wscale*w(:,i);
    end
  end
  
  % Create mdp copies of n
  nmtrx = repmat(n',1,mdp);

  iter_dir_max = 10;
  iter_mag_max = 10;
  tol = 1e-8;
  
  % Now begin the pursuit
  iter = 0;
  isdone = 0;
  X = x*w;
  [P,PCV,dPdX,dPdz] = pderivatives(X,n);
  while ~isdone
    bad_steps = 0;
    %
    % Calculate the length update
    %
    % Compute the figure-of-merit for overall magnitude
    lPcv = zeros(size(PCV));
    isOK = logical(zeros(size(PCV)));
    for idp = 1:mdp
      [lPcv(:,idp),isOK(:,idp)] = logpcvrobust(PCV(:,idp),...
                                               options.outlier_frac);
    end
    % The follow is -negentropy (the plusentropy?)
    % We want to minimize this.
    % Note it uses cross-validation
    Nvec = sum(isOK);
    Jmag = -(n*lPcv)./Nvec - log(diag(w'*w)')/2 - (1+log(2*pi))/2;
    PCV(~isOK) = 1;
    dznum = sum(nmtrx.*isOK.*dPdz./PCV)./Nvec + 1;
    dzdenom = sum(nmtrx.*isOK.*dPdz.^2./PCV.^2)./Nvec + 1;
    dz = dznum./dzdenom;
    for idp = 1:mdp
      is_smaller = 0;
      iter_mag = 0;
      step_size = 1;
      while ~is_smaller
        wtry = w(:,idp)*abs((1+step_size*dz(idp)));
        Xtry = x*wtry;
        [P(:,idp),PCVtry,dPdX(:,idp)] = pderivatives(Xtry,n);
        lPCVtry = log(PCVtry);
        lPCVtry(~isOK(:,idp)) = 0;  % any that were killed before should
                                    % be killed now
        Jmag_try = -(n*lPCVtry)/Nvec(idp) ...
          - log(wtry'*wtry)/2 - (1+log(2*pi))/2;
        is_smaller = (Jmag_try <= Jmag(idp) || step_size == 0);
        if ~is_smaller
          step_size = step_size/2;
          bad_steps = bad_steps+1;
        end
        iter_mag = iter_mag+1;
        if (iter_mag > iter_mag_max)
          step_size = 0;
          %warning('Failed to improve during the magnitude update')
        end
      end
      w(:,idp) = wtry;
    end


    %
    % Calculate the direction update
    %
    if options.oversmooth
      X = x*w/(1+options.oversmooth);
      [P,PCV,dPdX] = pderivatives(X,n);
    end
    % Compute the figure-of-merit for direction
    % This is also the negentropy, but here we are not using
    % cross-validation.
    Jdir = - (n*(isOK.*log(P)))./Nvec - log(diag(w'*w)')/2;
    % Compute the new direction
    slopelog = (dPdX./P).*isOK.*nmtrx;
    hesslog = sum((dPdX./P).*slopelog)./Nvec;
    for idp = 1:mdp
      % This next is -gradw(Jdir)
      gradwJdir(:,idp) = ...
          sum(repmat(slopelog(:,idp),1,options.maxcomponents).*x,1)/Nvec(idp) + ...
          w(:,idp)'/(w(:,idp)'*w(:,idp));
    end
    % Now adjust the direction to insure that the different directions in
    % w will be orthogonal (really Cx-conjugate, but we've made Cx = 1)
    wnorm2 = w'*w;   % w will not in general be a unit vector, because we
                     % use the magnitude to encode the amount of smoothing
    beta = wnorm2\(w'*gradwJdir);
    beta = (beta+beta')/2;  % symmetrize
    % Calculate the hessian. Since beta is the "current approximate
    % value," it's in fact possible for the naive hessian,
    %      diag(hesslog) + beta,
    % to not be positive definite.  One option is to simply throw
    % out the beta term.  The following is slower but (perhaps?) better?
    %hessJdir = diag(hesslog) + sqrtm(beta'*beta);
    hessJdir = spdiags(hesslog',0,mdp,mdp);  % better when dp is large
    dw = (gradwJdir - w*beta)/hessJdir;
    dwnorm2 = diag(dw'*dw);  % for convergence testing
    % OK, now we know what directions we're going to be looking in.
    % But here, we'll keep the vectors in w of the same length
    is_smaller = 0;
    step_size = 1;
    iter_dir = 0;
    while ~is_smaller
      if (iter+1 < length(options.random_magnitude))
        % We're doing stochastic descent
        % Add a random vector of size |dw|*options.random_magnitude(iter)
        % Note the length of dw is summed over all the directions, so for
        % directions with little gradient, the random step in that
        % direction can be comparatively big.  The thinking here is that
        % when the gradient isn't very steep, we can afford to move around.
        dw = dw + options.random_magnitude(iter+1)*...
             sqrt(mean(dw(:).^2))*randn(size(dw));
        wtry = w + dw;
        wtry = wtry/sqrtm(wtry'*wtry);  % Make them orthonormal
        wtry = wtry*diag(sqrt(diag(wnorm2)));  % Restore their overall length
        step_size = 1;
        is_smaller = 1;
        bad_steps = -options.random_magnitude(iter+1);
        if ~options.oversmooth
          % We have to prepare for the length optimization, but only if
          % we're not oversmoothing (because if so, we'll do it anyway down
          % below)
          Xtry = x*wtry;
          [Ptry,PCV,dPdXtry,dPdz] = pderivatives(Xtry,n);
        end
      else
        % We're greedy.  Update the direction in a guaranteed-convergent
        % way, making certain that the step taken decreases the
        % figure-of-merit.
        wtry = w + step_size * dw;
        wtry = wtry/sqrtm(wtry'*wtry);  % Make them orthonormal
        wtry = wtry*diag(sqrt(diag(wnorm2)));  % Restore their overall length
        Xtry = x*wtry/(1+options.oversmooth);
        [Ptry,PCV,dPdXtry,dPdz] = pderivatives(Xtry,n);
        Jdir_try = -(n*(isOK.*log(Ptry)))./Nvec - log(diag(wtry'*wtry)')/2;
        is_smaller = (sum(Jdir_try) <= sum(Jdir) || step_size == 0); % to handle roundoff errors
        if ~is_smaller
          step_size = step_size/10;
          bad_steps = bad_steps+1;
        end
        iter_dir = iter_dir+1;
        if (iter_dir > iter_dir_max)
          step_size = 0;
          %warning('Failed to improve during the direction update')
        end
      end
    end
    step_size_dir = step_size;
    w = wtry;
    % If we're oversmoothing, we have to leave the calculated values of
    % PCV, etc, in a non-oversmoothed state
    if options.oversmooth
      X = x*w;
      [P,PCV,dPdX,dPdz] = pderivatives(X,n);
    else
      X = Xtry;
    end
    
    switch options.display
      case 'full'
        disp(sort(Jmag))
        fprintf('.%g',bad_steps);
      case 'iter'
        fprintf('.');
    end
    
    % Are we done yet?
    iter = iter+1;
    dw2tot(iter) = sum((step_size_dir'.^2.*dwnorm2)');
    isdone = all(step_size_dir'.^2.*dwnorm2 < tol*diag(wnorm2)) || ...
             iter > options.itermax + length(options.random_magnitude);
  end
  if (iter > options.itermax + length(options.random_magnitude))
    warning('PPEK failed to converge');
  end
  if ~strcmp(options.display,'none')
    fprintf(' (%d iterations)\n',iter);
  end
  % Sort them in order of decreasing "goodness"
  [sJ,sort_order] = sort(Jmag); % since Jmag ~ -negentropy, want in increasing order
  w = w(:,sort_order);
  % Now invert the whitening to give a vector back in the correct direction
  w = (E*(D\w))';
  
  
  
function options = ppek_default_options(options,d)
  if ~isfield(options,'start')
    options.start = 'random';
  end
  if ~isfield(options,'dp')
    options.dp = d;
    if isnumeric(options.start)
      options.dp = size(options.start,1);
    end
  end
  options.dp = min(options.dp,d);
  if ~isfield(options,'maxcomponents')
    options.maxcomponents = d;
  end
  if ~isfield(options,'setmag')
    options.setmag = 0;
  end
  if ~isfield(options,'oversmooth')
    options.oversmooth = 0.5;
  end
  if ~isfield(options,'outlier_frac')
    options.outlier_frac = 0;
  end
  % The following is for stochastic descent
  if ~isfield(options,'random_magnitude')
    if (strcmp(options.start,'pca') || isnumeric(options.start))
      options.random_magnitude = [];
    else
      n_random = min(200,4*d);
      options.random_magnitude = 4*[ones(1,n_random/2) linspace(1,1/n_random,n_random/2)];
    end
  end
  if ~isfield(options,'itermax')
    options.itermax = min(1000,100+5*d);
  end
  if ~isfield(options,'display')
    options.display = 'iter';
  end
  