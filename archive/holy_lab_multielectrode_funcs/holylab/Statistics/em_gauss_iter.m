function [alphao,muo,sigma2o,p,nzIndex] = em_gauss_iter(X,n,alpha,mu,sigma2,options)
% EM_GAUSS_ITER: one iteration of EM algorithm
% Syntax:
%   [alpha,muout,sigma2out] = em_gauss_iter(X,n,alpha,mu,sigma2)
%   [alpha,muout,sigma2out,p] = em_gauss_iter(X,n,alpha,mu,sigma2)
%   [alpha,muout,sigma2out,p,nzIndex] = em_gauss_iter(X,n,alpha,mu,sigma2)
% where
%   X is a d-by-M matrix;
%   n is the multiplicity of each point of X (1-by-M);
%   alpha is a 1-by-q vector of gaussian amplitudes;
%   mu is a d-by-q matrix of gaussian centers (means);
%   sigma2 is a 1-by-q vector of gaussian variances;
% and
%   alphaout is a set of improved gaussian amplitudes;
%   muout is a set of improved gaussian centers;
%   sigma2out is a set of improved gaussian variances;
%   p is a N-by-q matrix of likelihoods, p(i,j) the likelihood that x_i
%     was drawn from the jth gaussian.  These p's are NOT normalized!
%   nzIndex is the index of the original components which came back with
%     non-zero amplitude;
%
% Also:
%   options is a structure which may have the following fields:
%     hard (default false): if true, does hard clustering
%       (p(i,j) = 0 or 1). If hard-clustering, the p output is instead an
%       index, p(i) is the most likely gaussian explaining the data point
%       x_i.  If a gaussian was discarded on this iteration, p will be
%       empty; you must iterate again to obtain a meaningful p.
%
% Note: when a particular gaussian explains less than two data points,
% the amplitude of that gaussian is sent to 0 and it is discarded. You
% can keep track of the discards, if necessary, using the nzIndex output.
%
% See also: EM_GAUSS.
  
% Copyright 2005 by Timothy E. Holy

  fast = 1;    % affects which algorithm is used for hard clustering

  if (nargin < 6)
    options = struct;
  end
  if ~isfield(options,'hard')
    options.hard = 0;
  end
  
  [d,M] = size(X);   % M = length(n), N = sum(n)
  N = sum(n);
  [d,q] = size(mu);
  % Calculate p, where
  %  p(j,i) = alpha(j)/(2*pi*sigma2(j))^(d/2) exp(-(x_i-mu_j)^2/2sigma2(j))
  sd = sqrdist(X,mu);   % This is N-by-q
  sigma2rep = repmat(sigma2,M,1);
  if (options.hard && fast)
    % Use log-likelihoods to avoid slow exponentiation
    lcoef = log(alpha) - d*log(sigma2)/2;
    ll = sd./sigma2rep;
    ll = ll - repmat(lcoef,M,1);  % This is -log-likelihood
    [minll,p] = min(ll,[],2);  % Minimize because it's -log-likelihood
    [clabel,n_per_gaussian] = agglabel(p);
    % Eliminate ones that account for fewer than 2 points
    isTooSmall = (n_per_gaussian < 1);
    clabel(isTooSmall) = [];
    n_per_gaussian(isTooSmall) = [];
    if any(isTooSmall), p = []; end   % Signal that p is not valid
    nzIndex = find(~isTooSmall);
    q = length(nzIndex);
    % Calculate the new amplitudes, means, and variances
    alphao = nan(1,q);
    muo = nan(d,q);
    sigma2o = nan(1,q);
    for i = 1:length(clabel)
      nclust = n(clabel{i});
      n_in_clust = sum(nclust);
      alphao(i) = n_in_clust/N;
      muo(:,i) = (X(:,clabel{i})*nclust')/n_in_clust;
      sigma2o(i) = (nclust*sd(clabel{i},nzIndex(i)))/n_in_clust;
    end
    return
  end
  coef = alpha./(2*pi*sigma2).^(d/2);
  sdsig = sd ./ sigma2rep;
  %p = repmat(coef,M,1).*exp(-sdsig/2); % p(i,j) refers to x_i & mu_j
  logp = repmat(log(coef),M,1) - sdsig/2;  % log-likelihood
  p = exp(logp);
  [logpmax,maxIndex] = max(logp,[],2);
  if options.hard
    error('This is messed up now')
    % Do hard clustering: each point explained by a single gaussian
    p = zeros(size(p));
    pindex = sub2ind(size(p),(1:M)',maxIndex);
    p(pindex) = exp(logpmax);
  end
  %pnorm = p ./ repmat(sum(p,2),1,q);  % fraction of x_i explained by j
  pnorm = exp(logp - repmat(logpmax,1,q));  % Robust version of pnorm calc.
  pnorm = pnorm ./ repmat(sum(pnorm,2),1,q);
  % Now we want to check to see if any components have become unused, as
  % evidenced by their amplitude becoming < 1/M
  alphatmp = sum(pnorm,1)/M;
  isTooSmall = (alphatmp < 2/M);  % these are ones we'll kill
  nzIndex = find(~isTooSmall);
  if any(isTooSmall)
    logp(:,isTooSmall) = [];   % Truncate any that explain less than 1 data point
    sd(:,isTooSmall) = [];
    p = exp(logp);
    q = size(p,2);
    % Now we have to re-normalize because we truncated some things
    [logpmax,maxIndex] = max(logp,[],2);
    pnorm = exp(logp - repmat(logpmax,1,q));  % Robust version of pnorm calc.
    pnorm = pnorm ./ repmat(sum(pnorm,2),1,q);
  end
  % Compute amplitude
  pnn = pnorm.*repmat(n',1,q);
  alphao = sum(pnn,1)/N; % note we introduced the multiplicities here
  % Compute the new means and variance
  muo = (X*pnn)./repmat(N*alphao,d,1);
  sigma2o = sum(pnn.*sd,1)./((d*N)*alphao);
