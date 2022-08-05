function [param2field,gfield2gparam,gf2gp_safe] = chainrule_gaussian(x,insertionFlag,normalizedFlag)
% CHAINRULE_GAUSSIAN: parametric Gaussians for optimization problems
%
% See CHAINRULE_WRAPPED_OPTFUNC for the general problem this helps solve.
% This defines a conversion between Gaussian parameters and a "field" (an
% array) of points.
%
% Syntax:
%   [param2field,gfield2gparam] = chainrule_gaussian(x,insertionFlag,normalizedFlag)
% where
%   x is a npoints-by-d array containing the coordinates of the positions
%     at which the Gaussian is to be evaluated
%   insertionFlag is an array of the size of A that acts as a
%     mask for the insertion of points in x into A.  If you want to use
%     all the points, just set equal to true(size(A)).
% and on output,
%   param2field is a function handle that converts Gaussian parameters to
%     the field
%   gfield2gparam is a function handle that converts gradients with
%     respect to the field into gradients with respect to the Gaussian
%     parameters.
%
% IMPORTANT NOTE: it is assumed that you will usually call gfield2gparam
% after the corresponding call to param2field.  This therefore uses
% nested functions (equivalent to persistent variables) to avoid
% recomputation of the gaussian field.  This, however, means that
% you will run into trouble if you want to call gfield2gparam for a
% different set of parameters than used in the most recent call to
% param2field.  In such cases, you should ask for a third output:
%   [param2field,gfield2gparam,gf2gp_safe] = chainrule_gaussian(x,insertionFlag)
%   gradp = gf2gp_safe(gradf,p)
% that will re-computate the field before applying the chain rule.
%
% Syntax of p:
%   p(1): the value of the gaussian at its peak (normalizedFlag = false),
%     or the integral in the continuum limit (normalizedFlag = true)
%   p(2:d+1): the centroid of the gaussian (in d-dimensional space)
% and one of the following:
%   p(d+2): 1/sigma^2, for an isotropic gaussian
%   p(d+2:2*d+1): 1./sigma.^2, for a diagonal covariance
%   p((1:d*(d+1)/2)+d+1): when the exponent of the gaussian is
%     (1/2)*x'*Omega*x, for a symmetric matrix Omega, this syntax
%     specifies each element of Omega. Convert between vector and
%     matrix forms of Omega using SQUAREFORM_DIAG.
%
% One interesting point to note is that you don't have to specify which
% of the forms for Omega you're planning to use. You can therefore do the
% optimization first with fewer parameters (e.g., a diagonal Omega) and
% then add the off-diagonal terms.
%
% A complete example: first, we define the geometry of the problem:
%   xr = 15; yr = 12;
%   R = 7; % points outside this circle are truncated to 0
%   coordArray = {-xr:xr,-yr:yr};
%   [X,Y] = ndgrid(coordArray{:});
%   insertionFlag = (X.^2 + Y.^2) < R^2;
%   x = [X(insertionFlag) Y(insertionFlag)];
% Now, create the function handles:
%   [param2field,gfield2gparam] = chainrule_gaussian(x,insertionFlag);
% Define a "target gaussian":
%   A0 = param2field([7 1 5 0.25]);
%   figure; imagesc(coordArray{:},insertionFlag)
%   figure; imagesc(coordArray{:},A0); colorbar
% Notice the gaussian is truncated by the circle. (You don't have to make
% it this complicated, this merely illustrates the possibilities.)
%
% Define an error function that works on arrays A:
%   function [val,grad] = errfunc(A,A0)
%     dA = A-A0;
%     val = sum(dA(:).^2)/2;
%     grad = dA;
% and set
%   func = @(A) errfunc(A,A0);
% (Notice this function is trivial, and would let you fit any functional
% form on any sized array.)
% 
% Define the full optimization function:
%   gfunc = @(p) chainrule_wrapped_optfunc(p,func,param2field,gfield2gparam);
% 
% Pick a (not very good) initial guess:
%    p0 = [1 0 0 1];
% Do the optimization:
%    p = conjgrad(gfunc,p0)
% and check the result.
%
% See also: CHAINRULE_WRAPPED_OPTFUNC, SQUAREFORM_DIAG.
  
% Copyright 2009 by Timothy E. Holy

  %% Input validation---do this once at initialization and then we don't
  %  have to worry about it anymore!
  sz = size(insertionFlag);
  sz1 = sz(sz > 1);  % Get rid of singleton dimensions
  if ~isnumeric(x)
    error('x must be numeric');
  end
  if (ndims(x) ~= 2)
    error(['x must be two-dimensional, with each row containing the' ...
	   ' coordinates of a single point']);
  end
  [npts,d] = size(x);
  if (length(sz1) ~= d)
    error('Size mismatch between x and insertionFlag');
  end
  if (nargin < 3)
    normalizedFlag = false;
  end
  if ~islogical(normalizedFlag)
    error('normalizedFlag must be true or false');
  end

  %% Define items that must be shared between the two nested functions and
  % thus act as "global" variables
  dx = [];
  dxOmega = [];
  Omega = [];
  diagOmega = [];
  g1 = [];
  g = [];
  psz = [];

  %% Define the function handles and return
  param2field = @(p) cg_param2field(p,x,insertionFlag,normalizedFlag);
  gfield2gparam = @(gradf) cg_gfield2gparam(gradf,insertionFlag,normalizedFlag);
  gf2gp_safe = @(gradf,p) cg_gf2gp_safe(gradf,p,x,insertionFlag,normalizedFlag);

  %% Nested functions
  function A = cg_param2field(p,x,insertionFlag,normalizedFlag)
    A = zeros(size(insertionFlag),class(x));
    d = size(x,2);
    psz = size(p);
    n_p = numel(p);
    amplitude = p(1);
    mu = p(2:d+1);
    normFactor = 1;
    if (n_p == d+2)
      Omega = p(d+2);
      if normalizedFlag
        normFactor = (Omega/(2*pi))^(d/2);
      end
    elseif (n_p == 2*d+1)
      diagOmega = p(d+2:end);
      Omega = spdiags(diagOmega(:),0,d,d);
      if normalizedFlag
        normFactor = sqrt(prod(diagOmega)/(2*pi)^d);
      end
    elseif (n_p == d*(d+1)/2+d+1)
      Omega = squareform_diag(p(d+2:end));
      if normalizedFlag
        normFactor = sqrt(abs(det(Omega))/(2*pi)^d);
      end
    else
      error('p does not have the correct size, given the dimensionality');
    end
    dx = x;
    for dimIndex = 1:d
      dx(:,dimIndex) = dx(:,dimIndex) - mu(dimIndex);
    end
    dxOmega = dx*Omega;
    dx2 = sum(dxOmega .* dx,2);
    g1 = exp(-dx2/2);  % the gaussian with unit amplitude
    if (normFactor ~= 1)
      g1 = normFactor*g1;
    end
    g = amplitude*g1;
    if (nargout > 0)
      A(insertionFlag) = g;
      A = reshape(A,size(insertionFlag));
    end
  end
    
  function gradp = cg_gfield2gparam(gradf,insertionFlag,normalizedFlag)
    gradp = zeros(psz);
    gradv = gradf(insertionFlag);
    % Gradient with respect to amplitude
    gradp(1) = sum(gradv .* g1);
    % Gradient with respect to mu
    for dimIndex = 1:d
      gradp(1+dimIndex) = sum(dxOmega(:,dimIndex) .* g .* gradv);
    end
    % Gradient with respect to parameters of Omega
    n_Omega = numel(gradp) - d - 1;
    if (n_Omega == 1 || n_Omega == d)
      gradtmp = zeros(1,d);
      if normalizedFlag
        if isscalar(Omega)
          diagOmega = Omega(ones(1,d));
        end
        for dimIndex = 1:d
          gradtmp(dimIndex) = sum((1/diagOmega(dimIndex)-dx(:,dimIndex).^2) .* g .* gradv);
        end
      else
        for dimIndex = 1:d
          gradtmp(dimIndex) = -sum(dx(:,dimIndex).^2 .* g .* gradv);
        end
      end
      if (n_Omega == 1)
        gradp(end) = sum(gradtmp)/2;
      else
        gradp(d+2:end) = gradtmp/2;
      end
    else
      % full-Omega version
      gradtmp = zeros(d,d);
      if normalizedFlag
        invOmega = inv(Omega);
        for dimIndex1 = 1:d
          for dimIndex2 = 1:dimIndex1
            gradtmp(dimIndex1,dimIndex2) = sum(g .* gradv .* ...
              (invOmega(dimIndex1,dimIndex2) - dx(:,dimIndex1) .* dx(:,dimIndex2)));
          end
        end
      else
        for dimIndex1 = 1:d
          for dimIndex2 = 1:dimIndex1
            gradtmp(dimIndex1,dimIndex2) = -sum(g .* dx(:,dimIndex1) .* ...
              dx(:,dimIndex2) .* gradv);
          end
        end
      end
      % Diagonals need adjustment because there aren't two entries per p
      for dimIndex = 1:d
        gradtmp(dimIndex,dimIndex) = gradtmp(dimIndex,dimIndex)/2;
      end
      gradp(d+2:end) = squareform_diag(gradtmp);
    end
  end
  
  function gradp = cg_gf2gp_safe(gradf,p,x,insertionFlag,normalizedFlag)
    cg_param2field(p,x,insertionFlag,normalizedFlag);
    gradp = cg_gfield2gparam(gradf,insertionFlag,normalizedFlag);
  end
end
