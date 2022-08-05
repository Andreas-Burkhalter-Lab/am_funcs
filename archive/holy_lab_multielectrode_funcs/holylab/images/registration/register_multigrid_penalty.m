function [imMg,detJnew,v,vgrad,hess] = register_multigrid_penalty(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,options)
% REGISTER_MULTIGRID_PENALTY: penalty, gradient, and hessian for nonrigid registration
% Syntax:
%   [imMg,detJnew,v,vgrad,hess] =
%   register_multigrid_penalty(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,options)
% where
%   u is the deformation (see REGISTER_MULTIGRID_VCYCLE for a full
%     explanation)
%   detJprev (default 1) is an array containing the determinant of a
%     "parent" deformation (useful for compositional improvement)
%   imM is the moving image
%   imF is the fixed image
%   mask is a logical array of the same size as imF, true in the pixels
%     for which you want to evaluate the data portion of the penalty
%   lambda is the coefficient for the log(det(J)) penalty (see
%     REGISTER_LOGDETPENALTY)
%   detJvalue is the target value for detJ (typically 1, unless you want
%     to scale the image)
%   upfunc is a function that prolongs the coarse-grid u to a fine-grid u
%     (see REGISTER_U_PROLONG; if u is supplied at full resolution, set to [])
%   downfunc is a function that restricts the full-resolution gradient
%     back down to the size of the supplied u (see
%     REGISTER_GRADU_RESTRICT; needed only if you are evaluating the
%     gradient)
%   options may have the following fields:
%     Mfull (default true) if true causes M(g(x)) to be evaluated in its
%       entirety; if false, it is evaluated only at points with mask = true
%     hessmax (default 10000): if the hessian has more rows than this, then
%       just evaluate the diagonal components.
%     zero_nans (default true): if true, any pixels that are NaN after
%       interpolation will be replaced by zero.
% and
%   imMg is the warped moving image
%   detJnew is the det(J) where J is the jacobian of u
%   v is the penalty value
%   vgrad is the penalty gradient with respect to u
%   hess is a nonnegative definite approximation to the hessian with
%     respect to u.
%
% See also: REGISTER_MULTIGRID_VCYCLE.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 10)
    options = struct;
  end
  options = default(options,'Mfull',true,'hessmax',10000,'zero_nans',true);
  if (nargout > 4)
    options.Mfull = true;
  end
  
  if (nargin < 9 || isempty(downfunc))
    downfunc = @(x) x;
  end
  if (nargin < 8)
    upfunc = @(x,y) copyinputs(x,y);
  end
  detJnew = 1;
  tmp = cell(1,nargout-1);
  if (lambda > 0)
    [ghdata,uhreg] = upfunc(u);
  else
    ghdata = upfunc(u);
  end
  [tmp{:}] = register_mismatch_noncov(ghdata,imM,imF,mask,options);
  imMg = tmp{1};
  p1 = tmp{2};
  if (nargout > 3)
    grad1 = tmp{3};
  end
  if (nargout > 4)
    dMdg = tmp{4};
  end
  if (lambda > 0)
    if (nargout <= 4)
      tmp = cell(1,nargout-1);
    else
      % We want the hessian. Use 5 outputs if diagonal-only, 6 if want the
      % full thing
      if (numel(uhreg) > options.hessmax)
        tmp = cell(1,5);
      else
        tmp = cell(1,6);
      end
    end
    if isscalar(detJprev)
      [tmp{:}] = register_logdetpenalty(uhreg,detJvalue/detJprev);
    else
      [tmp{:}] = register_logdetpenalty(uhreg,detJvalue,detJprev);
    end
    detJnew = tmp{1};
    p2 = tmp{2};
    if (nargout > 3)
      grad2 = tmp{3};
    end
    if (nargout > 4)
      if (numel(uhreg) > options.hessmax)
        % diagonal-only
        hess2 = sparse(tmp{4}+1,tmp{4}+1,tmp{5});
      else
        % full hessian
        hess2 = sparse(tmp{4}+1,tmp{5}+1,tmp{6});
      end
    end
  end
  v = p1;
  if (lambda > 0)
    v = v + lambda*p2;
  end
  if (nargout > 3)
    % Bring the grad back down to the size of u
    grad1 = downfunc(grad1);
    vgrad = grad1;
    if (lambda > 0)
      grad2 = downfunc(grad2);
      vgrad = vgrad + lambda*grad2;
    end
  end
  if (nargout > 4)
    % Calculate the Hessian
    dMdgH = downfunc(dMdg);
    N = numel(u);
    hess = (2/numel(imF))*spdiags(dMdgH(:).^2,0,N,N);
    if (lambda > 0)
      hess = hess + lambda*hess2;
    end
  end
  if isnan(v)
    %warning('Mismatch value is NaN; image may be shifted beyond field of view');
    v = Inf;
    if (nargout > 3)
      vgrad = zeros(size(vgrad));
    end
  end
  if (nargout > 3)
    if any(isnan(vgrad(:)))
      error('gradient should not have NaNs...');
    end
  end
end

function varargout = copyinputs(varargin)
  varargout = varargin;
end