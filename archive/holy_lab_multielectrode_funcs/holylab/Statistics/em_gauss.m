function [alpha,mu,sigma2,p,nzIndex] = em_gauss(X,varargin)
% EM_GAUSS: expectation-maximization (EM) algorithm with spherical gaussians
% Given a set of N data points X in d dimensions, model the probability
% density as a sum of spherically-symmetric gaussians with
% individually-optimized amplitude, center (mean), and width.
%
% Syntax:
%   [alpha,mu,sigma2] = em_gauss(X,q)
%   [alpha,mu,sigma2] = em_gauss(X,Y)
%   [alpha,mu,sigma2,p] = em_gauss(X,q or Y)
%   [...] = em_gauss(X,n,q or Y)
%   [...] = em_gauss(X,n,q or Y,options)
% where
%   X is a d-by-M matrix;
%   n is a 1-by-M vector, where n(i) is the multiplicity of x_i (if not
%     supplied, will be calculated by UNIQUE_DATA);
%   q is the initial number of gaussians used (final number may be fewer,
%     if some gaussians are driven to zero amplitude);
%   Y is a d-by-q matrix, giving the initial points for mu;
%   options is a structure which may have the following fields:
%     skipunique (default false): if true, will not call UNIQUE_DATA even
%       if n is not supplied (sets n = ones(1,M));
%     hard (default false): if true, performs hard-clustering (each point
%       assigned to a single gaussian rather than being explained by a
%       sum-of-gaussians);
%     movemu (default true): if true, permits the centers to move,
%       otherwise we just optimize amplitudes and widths;
%     ploteach (default false): if true, the data & centers are plotted
%       on each iteration;
%     plotpause (default false): if true, the algorithm pauses for
%       keyboard input after plotting;
%     plotproject: if supplied, projects the data & centers down to two
%       dimensions before plotting (a 2-by-d matrix) (default: use first
%       two dimensions of the input);
%     ptol (default 1/(10*sqrt(N))): the convergence criterion on the
%       estimated probability density, the relative change has to be
%       smaller than ptol*pdensity;
% and
%   alpha is a vector of amplitudes (1-by-ngaus);
%   mu is a matrix of centers (d-by-ngaus);
%   sigma2 is a vector of widths (1-by-ngaus);
%   p is a matrix (M-by-q), where p(i,j) is the likelihood that
%     x_i was drawn from the jth gaussian.  These p's are NOT normalized,
%     i.e., if you want to use this for soft clustering, you will first
%     want to normalize so that sum_j p(i,j) = 1.  However, the
%     un-normalized version might be useful for, e.g., computing the
%     density at the data points (pdensity = sum(p,2)).
%
% See also: KMEANS_HARD, EM_GAUSS_ITER.

% Copyright 2005 by Timothy E. Holy
  
  [d,N] = size(X);  % Note that N will be corrected later if nec.
  n = [];
  carg = 1;
  if (isvector(varargin{carg}) && length(varargin{carg}) == N)
    % X,n syntax
    n = varargin{1};
    carg = 2;
  end
  % This next might actually be Y, but em_gauss_init will figure it out
  % correctly
  q = varargin{carg};  
  if (length(varargin) > carg)
    options = varargin{carg+1};
  else
    options = struct;
  end

  options = em_gauss_options(options,d);

  % Reduce the input data to a unique set of points, with
  % multiplicities. We do this to facilitate the detection and
  % elimination of gaussians that are headed towards delta-functions
  % centered at a single position.
  if isempty(n)
    if options.skipunique
      n = ones(1,size(X,2));
    else
      [X,n] = unique_data(X');
      X = X';
    end
  end
  % Now correct N if necessary
  M = size(X,2);
  N = sum(n);

  if options.ploteach
    % Plot the data once; the centers will be plotted freshly on each
    % iteration
    xp = options.plotproject*X;
    hdata = plot(xp(1,:),xp(2,:),'.');
  end
  
  % Initialize the guesses
  [mu,sigma2,p,nzIndex] = em_gauss_init(X,q);
  q = size(mu,2);   % in case q was really Y, or points were tossed
  alpha = sum(p,1)/N;  % Note here that p = pnorm, so this is OK
  pdensity = sum(p,2);  % This will later be an estimate of the prob. density
  if ~isfield(options,'ptol')
    options.ptol = 1/(10*sqrt(N));
  end
  
  if options.ploteach
    hcircle = em_gauss_plotcenters(mu,sigma2,options);
  end
  
  % Iteratively improve the mixture of gaussians
  done = 0;
  iter = 0;
  while ~done
    pdensityold = pdensity;
    % Update the mixture
    [alpha,munew,sigma2,p,nzIndextmp] = em_gauss_iter(X,n,alpha,mu, ...
                                                      sigma2,options);
    nzIndex = nzIndex(nzIndextmp);
    if options.movemu
      mu = munew;
    elseif (length(nzIndextmp) < size(mu,2))
      mu = mu(:,nzIndextmp);
    end
    % Plot, if desired
    if options.ploteach
      hcircle = em_gauss_plotcenters(mu,sigma2,options,hcircle);
    end
    % Convergence criterion: that the density computed for the data set
    % has stopped changing
    % But first we have to deal with the weird case when using hard
    % clustering
    if isempty(p)
      iter = iter+1;
      continue
    end
    pdensity = sum(p,2);
    done = all(abs(pdensity-pdensityold) < options.ptol*pdensity);
    iter = iter+1;
    %pdensitychange(iter) = max(abs(pdensity-pdensityold)./pdensity);
  end
  %iter

function [mu,sigma2,p,keepIndex] = em_gauss_init(X,q)
  [d,M] = size(X);
  if isscalar(q)
    % q syntax: pick random starting points for mu
    rpIndex = randperm(M);
    mu = X(:,rpIndex(1:q));
  elseif (size(q,1) == d)
    % Y syntax: use the user's input for starting mu
    mu = q;
    q = size(mu,2);
  else
    error('Don''t recognize syntax for q/Y')
  end
  % Compute the other parameters by hard clustering on mu
  [dist,index] = mindist(X,mu);
  p = zeros(M,q);
  pindex = sub2ind(size(p),1:M,index);
  p(pindex) = 1;   % p = 1 for the closest center
  clabel = agglabel(index);
  sigma2 = zeros(1,q);
  for i = 1:q
    sigma2(i) = mean(dist(clabel{i})); % variance of pts belonging to center
  end
  keepIndex = find(sigma2);
  if any(sigma2 == 0)
    mu = mu(:,keepIndex);
    sigma2 = sigma2(keepIndex);
    [dist,index] = mindist(X,mu);
    p = zeros(M,length(keepIndex));
    pindex = sub2ind(size(p),1:M,index);
    p(pindex) = 1;   % p = 1 for the closest center
  end

  
function hcircle = em_gauss_plotcenters(mu,sigma2,options,hcircle)
  if (nargin > 3)
    delete(hcircle(ishandle(hcircle)))
  end
  % Plot the circles on each iteration
  mup = options.plotproject*mu;
  [d,q] = size(mu);
  hplot = sqrt(2*sigma2/d);  % estimate fraction of variance in coords 1&2
  theta = linspace(0,2*pi,200);
  ctheta = cos(theta);
  stheta = sin(theta);
  for i = 1:q
    hcircle(i) = line(mup(1,i)+hplot(i)*ctheta,...
                      mup(2,i)+hplot(i)*stheta,...
                      'Color','r');
    %hcenter = line(mup(1,:),mup(2,:),...
    %               'LineStyle','none',...
    %               'Marker','o',...
    %               'Color','r');
  end
  axis equal
  if options.plotpause
    pause
  else
    drawnow
  end

  
function options = em_gauss_options(options,d)
  if ~isfield(options,'ploteach')
    options.ploteach = 0;
  end
  if ~isfield(options,'plotpause')
    options.plotpause = 0;
  end
  if ~isfield(options,'plotproject')
    options.plotproject = [eye(2), zeros(2,d-2)];
  end
  if ~isfield(options,'skipunique')
    options.skipunique = 0;
  end
  if ~isfield(options,'hard')
    options.hard = 0;
  end
  if ~isfield(options,'movemu')
    options.movemu = 1;
  end
  % Note we don't do ptol here because we don't know N