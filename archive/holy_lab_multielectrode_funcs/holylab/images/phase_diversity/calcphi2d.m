function [phi,fval,object] = calcphi2d(image_series,H0,phi0,options)
% CALCPHI: calculate aberrations from phase-diverse images
%
% This function optimizes phase-aberrations to account for a set of
% images produced from a single (unknown) object.
%
% Syntax:
%   p = calcphi(image_series,H0,p0,options)
% where
%   image_series: an array of K 2d-images, h-by-w-by-K
%   H0: the pupil mask (h-by-w, or h-by-w-by-K if you need a different
%     pupil for each image)
%   p0: the initial guess for phi or the parameters that determine phi
%     (see below)
%   options: a structure which may have the following fields:
%     param2phi, gphi2gparam: functions which convert parameters to phi
%       and gradients in phi to gradients in the parameters,
%       respectively.  See CALCPHI_LINEARSUMS and CALCPHI_ZERNIKE.
%     normalize_grad (default false): if true, constrains the maximum
%       change in any parameter to be +/-mu_max (see below) on any
%       iteration cycle
%     mu_max: the largest linesearch step size accepted. If normalize_grad
%       is true, this defaults to 1.
%     Display (default true): print convergence information on each
%       iteration
%     use_smooothed_gradient (default false): if true, uses graddescent_smoothed
%       for optimization
%   You can also supply any other fields used by "conjgrad" (e.g., iter_max).
% and
%   p is the optimized set of parameters describing the phases. 
%
% If you don't supply the options structure, the only mode supported is
% when image 1 is unaberrated and image 2 is aberrated; in that case
% the full pupil phase will be optimized.
%
% See also: CALCPHI_LINEARSUMS, CALCPHI_ZERNIKE.
  
% Copyright 2009 by Timothy E. Holy
  
  %% Argument parsing
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'normalize_grad',false,'Display',true,'use_smoothed_gradient',false);
  if options.normalize_grad
    options = default(options,'mu_max',1);
  end
  
  impsz = size(image_series);
  if (length(impsz) ~= 3)
    error('Image data must be a 3-dimensional array')
  end
  imsz = impsz(1:2);
  K = impsz(3);
  M = size(H0,3);
  
  if (~isfield(options,'param2phi') || ~isfield(options,'gphi2gparam'))
    if (K ~= 2)
      error(['If you are supplying more than two images, you need to' ...
	     ' specify conversion functions param2phi and gphi2gparam']);
    end
    % With two images, we default to assuming that image1 is unaberrated
    % and image2 is aberrated.
    mixPhi = [0 1];
    s = calcphi_linearsums(mixPhi,0);
    options = copyfields(s,{'param2phi','gphi2gparam'},options);
  end

  % Initial guess for phi0
  if ((nargin < 3 || isempty(phi0)) && K == 2)
    % No initial guess for phi supplied, make it random and for a single
    % aberrated path
    phi0 = 0.1*randn(imsz).*H0(:,:,1);
  end
  
  % Create the function to be optimized
  optfunc = @(p) cp_pdwrapper(p,image_series,H0,options);
  
  % Do the optimization
  if options.use_smoothed_gradient
    minops = options;
    minops.mask = H0;
    schedule = [8 8; 4 4; 2 2; 1 1; 0 0];
    [phi,fval] = graddescent_smoothed(optfunc,phi0,schedule,options);
  else
    [phi,fval] = conjgrad(optfunc,phi0,options);
  end
  if (nargout > 2)
    [val,grad,object] = optfunc(phi);
  end
end

function [val,grad,obj] = cp_pdwrapper(p,image_series,H0,ops)
  impsz = size(image_series);
  K = impsz(3);
  phik = ops.param2phi(p);
  if (nargout > 1)
    if (nargout > 2)
      [val,gradphi,obj] = pdpenalty(phik,image_series,H0);
    else
      [val,gradphi] = pdpenalty(phik,image_series,H0);
    end
    grad = ops.gphi2gparam(gradphi);
    if ops.normalize_grad
      % Note! The next two lines mean that the derivative isn't right, but it
      % will make mu=1 be a reasonable size step (given that relevant
      % phases are of size unity)
      mxgrad = max(abs(grad(:)));
      grad = grad / mxgrad;
    end
  else
    val = pdpenalty(phik,image_series,H0);
  end
end
