function [w,h,err] = ssparse_nnmf1(M,options)
% SSPARSE_NNMF1: non-negative sparse matrix factorization with a single component
%
% Non-negative matrix factorization with a single component is particularly
% simple: all of the coordinates are independent, and so on each iteration
% it suffices to solve the linear equation and just set any coordinates
% with negative solutions to zero.
%
% This algorithm takes its input in the form of a "sparse matrix structure"
% parametrized in terms of the indices and values of the entries rather than a
% Matlab-format sparse matrix. (The Matlab matrices, for all their
% convenience for linear algebra, have some performance issues particularly
% with indexing and replacement.)
% The syntax is
%   [w,h] = ssparse_nnmf1(M)
%   [w,h,err] = ssparse_nnmf1(M,options)
% where
%   M is a "sparse matrix structure."
%
%   options may have the fields:
%     w0: an initial guess for w
%     h0: an initial guess for h (not used if you supply w0). This does not
%       have to be normalized
%     max_iter (default 10): the maximum number of alternating
%       optimizations
%     tol (default 1e-6): the tolerance for convergence
%
% On output, the matrix S corresponding to M is approximately represented
% as w*h. The convention is that h is normalized, sum(h.^2) = 1. err
% contains the summed square error, equal to sum(M(:).^2) - sum(w.^2).
%
% See also: NNMF.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'max_iter',10,'tol',1e-6);
  
  [index1,index2,value] = ssparse_find(M);
  index1 = index1(:);
  index2 = index2(:);
  value = value(:);
  maxi1 = max(index1);
  if isfield(options,'w0')
    w = options.w0(:);
  elseif isfield(options,'h0')
    w = accumarray(index1,options.h0(index2) .* value)/sqrt(sum(options.h0.^2));
  else
    w = ones(maxi1,1);
  end
  w(w<0) = 0;
  w2 = sum(w.^2);
  if (w2 == 0)
    w = ones(maxi1,1);
    w2 = maxi1;
  end
  
  w2Old = 0;
  iter = 0;
  while (1)
    h = accumarray(index2,w(index1) .* value) / w2;
    h(h < 0) = 0;
    h = h / sqrt(sum(h.^2));  % h is normalized
    w = accumarray(index1,h(index2) .* value);
    w(w < 0) = 0;
    w2 = sum(w.^2);
    iter = iter+1;
    if (w2 - w2Old < options.tol * w2 || iter > options.max_iter)
      break
    end
    w2Old = w2;
  end
  h = h';
  if (nargout > 2)
    err = sum(value.^2) - sum(w.^2);
  end
