function negen = negentropy(X,n,options)
% NEGENTROPY: calculate the negentropy of a projected data set
% Syntax:
%   negen = negentropy(X)
%   negen = negentropy(X,n)
%   negen = negentropy(X,n,options)
% where
%   X is a 1-by-M vector containing your (projected) data
%   n is a 1-by-M vector containing the multiplicities (defaults to ones)
% and
%   negen is the negentropy, H(p_gaussian) - H(p), where p is the kernel
%     estimate of the density (not MPL, so it's consistent with
%     projpurscg), and p_gaussian is a gaussian with the same variance as
%     the data.
% You can control the calculation by setting the following fields of
% options:
%   optim_length (default true): if true, calls mpl_optw to find the
%     optimimum projection factor (i.e., kernel width) for estimating the
%     density.  Set this to false only if X are already projected with the
%     optimum factor---note that the outputs of projpurscg are _not_
%     optimized (though they'll be close) because of oversmoothing and the
%     fact that usually the length optimization is terminated before the
%     direction optimization.
%   outlier_frac (default 0): helps with handling outliers, which can
%     distort the optimal W. Specifies the fraction of points to be
%     treated as outliers, and thus not contribute to the
%     cross-validated action. 0 retains all points; 0.05 would discard
%     the 5% least-likely points.
%
% See also: PPEK, PROJPURSCG.
  
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'optim_length')
    options.optim_length = 1;
  end
  
  if (size(X,1) > 1)
    error('negentropy works only in d=1');
  end
  if (nargin < 2)
    [X,n] = uniquemult(X');
    X = X';
  end
  [d,M] = size(X);
  if (length(n) ~= M)
    error('The sizes of X and n do not match');
  end
  N = sum(n);
  
  ops = struct('kernel_only',1,'display','none');
  if isfield(options,'outlier_frac')
    ops.outlier_frac = options.outlier_frac;
  end
  if options.optim_length
    p = mpl_optw(X,n,ops);
  else
    p = mpl(X,n,1,ops);
  end
  H = -(n*log(p))/N;
  
  % Compute covariance matrix of input points, properly accounting for
  % multiplicity
  nd = spdiags(n',0,M,M);
  Xm = mean(nd*X');
  dX = X - Xm;
  CX = (dX*nd*dX')/(N-1);
  Hg = (1/2)*log(CX) + (1 + log(2*pi))/2;

  negen = Hg - H;
  