function x = irls(A,b,x,frac)
% IRLS: iteratively-rescaled least squares (solves least abs. dev)
% This function solves the least absolute deviation problem, a robust
% version of least squares.  The problem is
%     min |A*x - b|
% where | | is the L1-norm.  IRLS solves this problem by rescaling the
% least-squares equation: if r is the vector of residuals, i.e. r = b -
% A*x, then IRLS defines a new problem
%    min (W*A*x - W*b)^2
% where the re-weighting matrix W is approximately W = diag(1./|r|).
% Syntax:
%   xout = irls(A,b)
% This will give you the solution starting from the least-squares guess.
%
% If you want finer-grained control,
%   xout = irls(A,b,x,frac)
% where
%    x is your initial guess for the solution (omit or supply [] if you
%      want to calculate it by least-squares),
%    frac is a vector describing a "regularization sequence" for
%      solving the problem:  on each iteration, residuals |r| that are of
%      smaller magnitude than a threshold equal to frac(i)*max(|r|) are
%      set to the threshold for the purpose of W.  This prevents division
%      by zero or other similar problems.  The exact solution is obtained
%      in the limit frac -> 0; you supply a vector of values to try.
%      Default is frac = [1e-1 1e-2 1e-3 1e-4].

% Copyright 2007 by Timothy E. Holy
  
  if (nargin < 3)
    x = [];
  end
  if (nargin < 4)
    frac = [0.1 0.01 0.001 0.0001];
  end
  if isempty(x)
    x = A\b;
  end
  n_iter = length(frac);
  for i = 1:n_iter
    r = A*x-b;
    w = abs(r);
    thresh = frac(i)*max(w);
    w(w < thresh) = thresh;
    w = thresh./w;
    WA = diag(w)*A;
    Wb = b .* w;
    x = WA\Wb;
  end
