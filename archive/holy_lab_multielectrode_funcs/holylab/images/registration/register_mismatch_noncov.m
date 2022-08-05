function [imMg,val,grad,dMdg] = register_mismatch_noncov(g,imM,imF,mask,options)
% REGISTER_MISMATCH_NONCOV: compute the mismatch error under deformation
%
% This function implements warping and penalty-function evaluation for
% non-covariant image registration. The penalty is
%      E = sum(mask(x) .* dI(x).^2)
% where dI = imFixed - imMoving(g(x)) and mask is a logical array of the
% same size as the images, 1 in places where you want to include the mismatch.
%
% Syntax:
%   [imMg,val] = register_mismatch_noncov(g,imM,imF,mask)
%   [imMg,val] = register_mismatch_noncov(g,imM,imF,mask,options)
% This calculates the deformed moving image imMg, evaluating it at the
% positions defined by the array g (with size [size(imMg) ndims(imM)]).  It
% also calculates penalty, supplying the fixed image imF and the mask.
%   The two syntaxes differ in an important way: by default, imMg is
% calculated only at positions for which mask == true.  This can save
% considerable computation time.  If, however, you want to calculate the
% "full" imMg, then use the second syntax with Mfull = true.
%
%   [imMg,val,grad] = register_mismatch_noncov(g,imM,imF,mask)
% This version computes the gradient of the penalty with respect to g.
%
%   [imMg,val,grad,dMdg] = register_mismatch_noncov(g,imM,imF,mask)
% This version also returns the gradient of the moving image itself with
% respect to g. This can be valuable in approximating the Hessian (and does
% not require any extra computation).

% Copyright 2009-2010 by Timothy E. Holy
  
  if (nargin < 5)
    options = struct;
  end
  options = default(options,'Mfull',false,'zero_nans',true);

  szM = size(imM);
  dimFlag = szM > 1;
  if ~dimFlag(1)
    error('First dimension must not be unity (try supplying the transpose image)');
  end
  n_dims = sum(dimFlag);
  szg = size(g);
  szg(end+1:n_dims+1) = 1;
  if (szg(n_dims+1) ~= n_dims)
    error('Dimensionality of g does not match image');
  end
  
  if (nargout < 3)
    out = cell(1,1);
  else
    out = cell(1,2);
  end
  
  N = sum(mask(:));
  if options.Mfull
    % Compute the full-size interpolated image
    [out{:}] = imqinterp(g,imM);
    imMg = out{1};
    if options.zero_nans
      imMg(isnan(imMg)) = 0;
    end
    dI = imF - imMg;
    val = sum(dI(mask(:)).^2)/N;
    if (nargout > 2)
      dMdg = out{2};
      notmask = ~mask;
      dI(notmask) = 0;
      killFlag = repmat(notmask,[ones(1,n_dims) n_dims]);
      dMdg(killFlag) = 0;  % kill NaNs in regions outside the mask
      grad = repmat((-2/N)*dI,[ones(1,n_dims) n_dims]) .* dMdg;
    end
  else
    % Compute the interpolated image just at the pixels with mask == true
    % (this can sometimes reduce computation time)
    n_elem = prod(szg(1:end-1));
    grs = reshape(g,[n_elem n_dims]);
    gmask = grs(mask(:),:);
    [out{:}] = imqinterp(gmask,imM);
    imMg = out{1};
    if options.zero_nans
      imMg(isnan(imMg)) = 0;
    end
    imF = imF(mask(:));
    dI = imF - imMg;
    val = sum(dI.^2)/N;
    if (nargout > 2)
      dMdg = out{2};
      gradLine = repmat((-2/N)*dI,[1 n_dims]) .* dMdg;
      % Must put grad into shape of g, filling in with zeros in places
      % where mask == false
      grad = zeros(szg);
      tmp = zeros(szg(1:end-1));
      colons = repmat({':'},1,n_dims);
      for dimIndex = 1:n_dims
        tmp(mask(:)) = gradLine(:,dimIndex);
        grad(colons{:},dimIndex) = tmp;
      end
    end
  end
