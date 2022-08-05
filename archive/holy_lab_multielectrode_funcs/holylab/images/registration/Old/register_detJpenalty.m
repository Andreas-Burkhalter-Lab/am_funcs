function [detJnew,val,grad] = register_detJpenalty(g,c,detJprev)
% REGISTER_DETJPENALTY: penalty for deviations from constant volume
%
% Given a deformation g, this function calculates the value and gradient of
% a penalty
%     E = sum_x (detJ - c)^2/N
% where the sum_x is over all pixels x, and J(x) is the jacobian matrix
% assciated with the deformation g(x). c is the target value for the
% determinant (c = 1 if you don't need to have an overall scale). N is the
% total # of points in the grid.
%
% Syntax:
%   detJ = register_detJpenalty(g)
% This version just computes detJ
%
%   [detJ,val] = register_detJpenalty(g,c)
% This also computes the value of the penalty function.
%
%   [detJ,val,grad] = register_detJpenalty(g,c)
% This also computes the gradient of the penalty function.
%
%   [detJnew,val,grad] = register_detJpenalty(g,c,detJprev)
% This version is particularly useful for multigrid, and implements the
% modified penalty
%     E = sum_x(detJprev.*detJnew - c)^2
% where only detJnew depends on g. detJprev might be the "cumulative
% jacobian" from all finer grids.
%
% This needs to be implemented as a MEX file if it is to have reasonable
% performance.

% Copyright 2009 by Timothy E. Holy

  if (isempty(g) && nargin < 3)
    if (c ~= 1)
      error('Empty input not supported for c other than c = 1');
    end
    detJnew = 1;
    val = 0;
    grad = 0;  % not a cell array, what to do?
    return
  end

  if ~isempty(g)
    n_dims = length(g);
    sz = size(g{1});
    npix = prod(sz);
    J = register_jacobian(g);
    detJnew = zeros(sz,class(g{1}));
    if (nargin < 3)
      detJprev = 1;%ones(sz,class(g{1}));
    end
    
    for i = 1:npix
      detJnew(i) = det(J(:,:,i));
    end
    if (nargout < 2)
      return
    end
    Jprod = detJprev.*detJnew;
    if any(Jprod(:) * c < 0)
      val = inf;
      return
    end
    dJ = Jprod - c;
    N = numel(dJ);
    val = sum(dJ(:).^2)/N;
    if (nargout < 3)
      return
    end
    M = zeros(size(J),class(J));
    coef = dJ .* Jprod;
    for i = 1:npix
      Ji = inv(J(:,:,i));
      M(:,:,i) = coef(i) * Ji;
    end
  else
    % Empty g but non-empty detJprev
    detJnew = detJprev;
    N = numel(detJprev);
    val = sum((detJprev(:) - c).^2);
    n_dims = ndims(detJprev);
    sz = size(detJprev);
    Mbase = (detJprev - c).*detJprev;
    M = zeros([n_dims,n_dims,size(detJprev)],class(detJprev));
    for i = 1:n_dims
      M(i,i,:) = Mbase(:);
    end
  end
  grad = cell(1,n_dims);
  for i = 1:n_dims
    grad{i} = zeros(sz,class(detJprev));
  end
  colons = repmat({':'},1,n_dims);
  for i = 1:n_dims
    for k = 1:n_dims
      if (sz(k) == 1)
        continue
      end
      Mtmp = reshape(M(k,i,colons{:}),size(detJprev));
      dM = -2*deriv(Mtmp,k);
      % Now fix the edges
      cc = colons;
      ccp = colons;
      ccm = colons;
      indxv = [1 2 sz(k)-1 sz(k)];
      indxv = indxv(indxv > 0 & indxv <= sz(k));
      for indx = indxv
        cc{k} = indx;
        ccp{k} = indx+1;
        ccm{k} = indx-1;
        tmp = zeros(size(Mtmp(cc{:})));
        if (indx > 1)
          tmp = Mtmp(ccm{:}) * (2 - (indx > 2));
        end
        if (indx < sz(k))
          tmp = tmp - Mtmp(ccp{:}) * (2 - (indx < sz(k)-1));
        end
        dM(cc{:}) = tmp - 2 * Mtmp(cc{:}) * ((indx < sz(k)) - (indx > 1));
      end
      grad{i} = grad{i} + dM;
    end
  end
  for i = 1:n_dims
    grad{i} = grad{i}/N;
  end
  