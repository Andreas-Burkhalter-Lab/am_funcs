function varargout = em_relax(X,weights,centroids,C,varargin)
% em_relax: optimize gaussian mixture model by EM, given a starting guess
%
% Syntax:
%   [weightsOut,centroidsOut,Cout] = em_relax(X)
% (syntax when you want a single gaussian for the whole data set)
%   [weightsOut,centroidsOut,Cout] = em_relax(X,weightsIn,centroidsIn,Cin)
%   [...,g,ll] = em_relax(X,...)
%   [...] = em_relax(...,options)
% where
%   X is the initial d-by-N data matrix, each d-dimensional point stored in
%     one column of X
%   weightsIn is a 1-by-k vector, containing the amplitude of each gaussian
%     in the initial mixture (non-negative, will be normalized so overall
%     scale is irrelevant)
%   centroidsIn is the d-by-k matrix of intial centroid positions
%   Cin is a d-by-d-by-k array, containing the initial guesses for the
%     covariance matrices of the gaussians
% and
%   weightsOut is 1-by-k and contains the optimized amplitude (sums to 1)
%   centroidsOut is d-by-k and contains the final matrix of centroid
%     positions
%   Cout is d-by-d-by-k and contains the covariance matrices; note these
%     are "regularized", see covariance_regularized.
%   g is k-by-N, containing the value of each gaussian evaluated at each
%     data point.
%   ll is the total log-likelihood (a scalar)
%
% The algorithm terminates after the "E" step, and the gaussian parameters
% are evaluated on the "M" step, so there may be small discrepancies
% compared to a direct evaluation. If this is problematic, simply decrease
% the tolerance using options, which has the following fields:
%   iter_max (default 100): the maximum number of iterations
%   tol (default 1e-4*N): tol is the absolute change in log-likelihood
%     required for convergence
%   meanCdiag_coef (default 1e-6): add this fraction of the mean diagonal
%     of C to each gaussian's C; this prevents complete singularity blow-ups.
%
% See also: kmeans_relax, covariance_regularized.

% Copyright 2010 by Timothy E. Holy

  %% Initialization
  options = struct;
  curarg = 1;
  while (length(varargin) >= curarg)
    if islogical(varargin{curarg})
      error('Unsupported syntax');
%       outlierFlag = varargin{curarg};
%       X = X(:,outlierFlag);
%       haveOutlierFlag = true;
    elseif isstruct(varargin{curarg})
      options = varargin{curarg};
    end
    curarg = curarg+1;
  end
  [d,N] = size(X);
  if (nargin < 3 || isempty(centroids))
    centroids = zeros(d,0);
  end
  n_g = size(centroids,2);
  if (n_g < 2)
    % With a single gaussian, we know the solution
    varargout = cell(1,max(3,nargout));
    varargout{1} = 1;
    centroids = mean(X,2);
    dX = bsxfun(@minus,X,centroids);
    [~,C] = covariance_regularized(dX);
    varargout{2} = centroids;
    varargout{3} = C;
    if (nargout > 3)
      d2 = sum(dX.*(C\dX),1);
      g = exp(-d2/2)/((2*pi)^(d/2)*sqrt(det(C)));
      varargout{4} = g;
      if (nargout > 4)
        varargout{5} = sum(log(g));
      end
    end
    return
  end
  
  if (nargin < 5)
    options = struct;
  end
  options = default(options,'tol',1e-4*N,'iter_max',100,'meanCdiag_coef',1e-6);
  g = zeros(n_g,N);
  llOld = -inf;
  
  for iter = 1:options.iter_max
    %% E step
    % Calculate the value of all gaussians at all data points
    for i = 1:n_g
      dX = bsxfun(@minus,X,centroids(:,i));
      d2 = sum(dX.*(C(:,:,i)\dX),1);
      g(i,:) = weights(i)*exp(-d2/2)/((2*pi)^(d/2)*sqrt(det(C(:,:,i))));
    end
    % Compute the normalized contribution. We have to be particularly
    % careful about outlying data points for which the gaussians are all
    % zero. This generally should not happen, unless the user supplies a
    % poor starting guess.
    gsum = sum(g,1);
    isnz = gsum > 0;
    p = zeros(size(g));
    p(:,isnz) = bsxfun(@rdivide,g(:,isnz),gsum(isnz));
    gsum(~isnz) = min(gsum(isnz));

    %% Test for convergence
    ll = sum(log(gsum));
    if (ll - llOld < options.tol)
      break
    end
    llOld = ll;
  
    %% M step
    weights = sum(p,2)/N;
    meanCdiag = zeros(d,1);
    for i = 1:n_g
      if (weights(i) == 0)
        continue
      end
      centroids(:,i) = sum(bsxfun(@times,p(i,:),X),2)/(N*weights(i));
      dX = bsxfun(@minus,X,centroids(:,i));
      [~,C(:,:,i)] = covariance_regularized(dX,p(i,:));
      meanCdiag = meanCdiag + diag(C(:,:,i));
    end
    % Protect against blow-ups by adding a small component of the mean
    % diagonal of C to all gaussians
    meanCdiag = diag(meanCdiag/n_g);
    C = bsxfun(@plus,C,options.meanCdiag_coef*meanCdiag);
  end
  
  %% Remove any "empty" gaussians
  keepFlag = weights >= 1/N;
  varargout = cell(1,max(nargout,3));
  varargout{1}= weights(keepFlag);
  varargout{2} = centroids(:,keepFlag);
  varargout{3} = C(:,:,keepFlag);
  if (nargout > 3)
    varargout{4} = g(keepFlag,:);
  end
  if (nargout > 4)
    varargout{5} = ll;
  end
