function [W,H,convergence_info] = nonneg_matrix_factorization(A,n_factors,options)
% NONNEG_MATRIX_FACTORIZATION: write A approximately as WH
% A matrix A may be decomposed as
%       A = WH
% where the decomposition is approximate, in a least-squares sense. Either
% W or H, or both, may be constrained to be nonnegative.
% Syntax:
%   [W,H] = nonneg_matrix_factorization(A,n_factors)
%   [W,H] = nonneg_matrix_factorization(A,n_factors,options)
% where
%   A is the data matrix
%   n_factors is the number of factors (# columns in W/rows in H)
%   options may have the following fields:
%     Wnn (default true): if true, W is constrained to be nonnegative
%     Hnn (default true):     "  , H            "
%     Winit: an initial value for W
%     Wrand (default false): prefer a random start for W, rather than SVD,
%       when Wnn is false (will always use random start if Wnn is true, unless Winit is supplied)
%     Hsparse: coefficient in front of |H|1^2 to make H sparse
%     tol (default 1e-5): the criterion for convergence
%     maxiter (default 50): the maximum number of iterations
%     display (default true): if true, causes progress to be printed

% Copyright 2007 by Timothy E. Holy

  [m,n] = size(A);
  
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'Wnn',true);
  options = default(options,'Hnn',true);
  options = default(options,'Hsparse',0);
  options = default(options,'tol',1e-5);
  options = default(options,'maxiter',50);
  options = default(options,'display',true);
  options = default(options,'Wrand',false);
  
  if isfield(options,'Winit')
    W = options.Winit;
    if (size(W,2) < n_factors)
      W = [W rand(m,n_factors-size(W,2))];
    end
  elseif (~options.Wnn && ~options.Wrand)
    [U,S,V] = svd(A,'econ');
    for i = 1:n_factors
      W(:,i) = U(:,i) * sign(median(V(:,i)));
    end
  else
    W = rand(m,n_factors);
    Wsum = sum(W);
    W = W ./ repmat(Wsum,[size(W,1) 1]);
  end
  %[Atmp,Wtmp] = sparsify(A,W,options.Hsparse);
  %H = nnmf_solve(Wtmp,Atmp,options.Hnn);
  % Don't apply sparseness yet
  H = nnmf_solve(W,A,options.Hnn);
  err = inf;
  iter = 0;
  store_info = false;
  if (nargout > 2)
    store_info = true;
    convergence_info.err = [];
    convergence_info.converged = true;
  end
  while (1)
    R = W*H - A;  % the residual (not to be confused with the Cholesky decomp)
    errold = err;
    R2 = R.^2;
    err = sum(R2(:)) + options.Hsparse*sum(sum(H).^2);
    if options.display
      fprintf('  NMF: iter %d, err = %g\n',iter,err);
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
    % Solve for W
    W = nnmf_solve(H',A',options.Wnn)';
    % Normalize W
    if options.Wnn
      Wsum = sum(W);
    else
      %[U,S,V] = svd(W,'econ');
      %W = U;
      Wsum = sqrt(sum(W.^2));
    end
    W = W ./ repmat(Wsum,[size(W,1) 1]);
    % Solve for H
    [Atmp,Wtmp] = sparsify(A,W,options.Hsparse);
    H = nnmf_solve(Wtmp,Atmp,options.Hnn);
    iter = iter+1;
  end
  if (iter > options.maxiter)
    if store_info
      convergence_info.converged = false;
    else
      warning('Failed to converge');
    end
  end
  
function X = nnmf_solve(M,A,isnn)
  n_cols = size(A,2);
  X = nan(size(M,2),n_cols);
  if isnn
    % Nonnegative solution
    % Perhaps this could be improved by pre-computing all the possible
    % Cholesky decompositions, as long as the combinatorics work out
    % favorably (e.g., very few factors).
    use_fnnls = false;
    %tic
    if use_fnnls
      X = fnnls(M'*M,M'*A);
    else
      for i = 1:n_cols
        X(:,i) = snnls(M,A(:,i));
      end
    end
    %toc
  else
    % Unconstrained solution
    [R,p] = chol(M'*M);
    if (p == 0)
      X = R\(R'\(M'*A));
    else
      p
      X = pinv(M'*M)*(M'*A);
    end
  end

function [A,W] = sparsify(A,W,sparsecoef)
if sparsecoef
  A = [A; zeros(1,size(A,2))];
  W = [W; sqrt(sparsecoef)*ones(1,size(W,2))];
end