function [T,W,convergence_info] = robust_nmf(X,n_factors,options)
% ROBUST_NMF: write X approximately as T*W
% A matrix X may be decomposed as
%       X = T*W
% where the decomposition is approximate, in a robust sense. 
% Data points are on columns of X. T = templates, W = weights
  
% The following needs to be updated.
% Either
% W or T, or both, may be constrained to be nonnegative.
% Syntax:
%   [W,T] = nonneg_matrix_factorization(X,n_factors)
%   [W,T] = nonneg_matrix_factorization(X,n_factors,options)
% where
%   X is the data matrix
%   n_factors is the number of factors (# columns in W/rows in T)
%   options may have the following fields:
%     Wnn (default true): if true, W is constrained to be nonnegative
%     Tnn (default true):     "  , T            "
%     Tinit: an initial value for W
%     Trand (default false): prefer a random start for T, rather than SVD,
%       when Tnn is false (will always use random start if Tnn is true, unless Winit is supplied)
%     Tsparse: coefficient in front of |T|1^2 to make T sparse
%     tol (default 1e-5): the criterion for convergence
%     maxiter (default 50): the maximum number of iterations
%     display (default true): if true, causes progress to be printed

% Copyright 2007 by Timothy E. Holy

  [m,n] = size(X);
  
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'Wnn',true);
  options = default(options,'Tnn',false);
  options = default(options,'Tsparse',0);
  options = default(options,'tol',1e-5);
  options = default(options,'maxiter',50);
  options = default(options,'display',true);
  options = default(options,'Trand',false);
  
  if options.Tnn
    error('Not implemented')
  end

  if isfield(options,'Tinit')
    T = options.Tinit;
    if (size(T,2) < n_factors)
      T = [T rand(m,n_factors-size(T,2))];
    end
  elseif (~options.Tnn && ~options.Trand)
    [U,S,V] = svd(X,'econ');
    for i = 1:n_factors
      T(:,i) = U(:,i) * sign(median(V(:,i)));
    end
  else
    T = rand(m,n_factors);
    Tsum = sum(T);
    T = T ./ repmat(Tsum,[size(T,1) 1]);
  end

  W = rnmf_solveW(T,X,options.Wnn);
  err = inf;
  iter = 0;
  store_info = false;
  if (nargout > 2)
    store_info = true;
    convergence_info.err = [];
    convergence_info.converged = true;
  end
  while (1)
    R = T*W - X;  % the residual (not to be confused with the Cholesky decomp)
    errold = err;
    err = sum(sqrt(sum(R.^2,1)),2);
    if options.display
      fprintf('  RNMF: iter %d, err = %g\n',iter,err);
    end
    if store_info
      convergence_info.err(end+1) = err;
    end
    if (abs(errold - err) < options.tol*(err + errold))
      break
    end
    if (iter > options.maxiter)
      break
    end
    % Solve for T
    T = rnmf_solveT(W,X);
    % Normalize T
    if options.Tnn
      Tsum = sum(T);
    else
      Tsum = sqrt(sum(T.^2));
    end
    T = T ./ repmat(Tsum,[size(T,1) 1]);
    % Solve for W
    W = rnmf_solveW(T,X,options.Wnn);
    iter = iter+1;
  end
  if (iter > options.maxiter)
    if store_info
      convergence_info.converged = false;
    else
      warning('Failed to converge');
    end
  end
  
function W = rnmf_solveW(T,X,isnn)
  n_cols = size(X,2);
  if isnn
    % Nonnegative solution
    % Perhaps this could be improved by pre-computing all the possible
    % Cholesky decompositions, as long as the combinatorics work out
    % favorably (e.g., very few factors).
    W = nan(size(T,2),n_cols);
    for i = 1:n_cols
      W(:,i) = snnls(T,X(:,i));
    end
  else
    % Unconstrained solution
    R = chol(T'*T);
    W = R\(R'\(T'*X));
  end

function T = rnmf_solveT(W,X)
  W = W';
  X = X';
  n_dims = size(X,2);
  T = nan(size(W,2),n_dims,class(X));
  for i = 1:n_dims
    T(:,i) = irls(W,X(:,i));
  end
  T = T';
  