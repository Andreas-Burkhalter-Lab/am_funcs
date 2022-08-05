function [Zcoefs,fval] = register_zernike(im,pupildata,params)
% register_zernike: shift (and optionally adjust higher-order aberrations) to "align" images
% 
% This function performs rigid image registration using fourier optics and
% a Zernike basis. Alternatively, you can use this to fit Zernike
% polynomials to try to describe a whole image series, e.g., if you're moving a
% single actuator, in terms of aberrations of a "base" image.
%
% Syntax:
%   [Zcoefs,fval] = register_zernike(im,pupildata,params)
%   [Zcoefs,fval] = register_zernike(im,params)
% where
%   im is a h-by-w-by-K series of images (for "self-registration"), or
%     alternatively a 1-by-2 cell array of such image series (for
%     "cross-registration")
%   pupildata is a structure which must have the following fields:
%     H0,rho,theta: pupil parameters (mask, radius, angular)
%   (in the second syntax, these fields are supplied as part of params)
%   params may have the following fields:
%     Zindex (default [1 2]): a list of Zernike coefficients
%       (single-indexing scheme)
%     Zc0 (default [0 0], or [0 0 randn(1,n)] for more than 2
%       coefficients): a starting guess for the Zernike coefficients
%     baseIndex (default 1): the index of the frame that is to be used as
%       the reference image
% and
%   Zcoefs is a 2-by-K (or 3-by-K) matrix giving the Zernike coefficients
%     for maximum alignment (in the case of K = 2 this will be converted
%     into a vector describing just the "test" image)
%   fval is a vector of mismatches (will be converted to a scalar for K=2)
%
% See also: REGISTER_FOCUS_ZERNIKE.

% Copyright 2009 by Timothy E. Holy

  %% Argument parsing
  register_mode = 'self';
  if iscell(im)
    register_mode = 'cross';
    imsz = size(im{1});
  else
    imsz = size(im);
  end
  imsz(end+1:3) = 1;
  K = imsz(3);
  if (K == 1 && strcmp(register_mode,'self'))
    error('Must supply multiple images');
  end
  if (nargin == 2)
    params = pupildata;
  else
    % Copy the pupildata over to params
    params.H0 = pupildata.H0;
    params.rho = pupildata.rho;
    params.theta = pupildata.theta;
  end
  params = default(params,'Zindex',[1 2],'baseIndex',1);
  if ~isfield(params,'Zc0')
    params.Zc0 = [0 0 randn(1,length(params.Zindex)-2)];
    %params.Zc0 = zeros(1,length(params.Zindex));
  end
    
  %% Configuration for pairwise comparisons
  params.mode = 'Zn pair';
  zf = phi_parametrizations(params);
  zf.normalize_grad = false;
  zf.tol = 1e-8;
  
  %% Do the optimization; compare each to the base image, using the
  %  previous result to initialize the next computation
  Zcoefs = repmat(params.Zc0,K,1);
  fval = zeros(1,K);
  % Moving to larger frame #s
  start = params.Zc0;
  indxStart = params.baseIndex;
  if strcmp(register_mode,'self')
    indxStart = params.baseIndex+1;
  end
  for indx = indxStart:K
    rz_work(indx);
  end
  % Moving to smaller frame #s
  if strcmp(register_mode,'self')
    Zcoefs(params.baseIndex,:) = 0;
    start = -Zcoefs(params.baseIndex+1,:);
  else
    start = Zcoefs(params.baseIndex,:);
  end
  for indx = params.baseIndex-1:-1:1
    rz_work(indx);
  end
  
  function rz_work(indx)
    switch register_mode
      case 'self'
        im_pair = im(:,:,[params.baseIndex indx]);
      case 'cross'
        im_pair = im{1}(:,:,indx);
        im_pair(:,:,2) = im{2}(:,:,indx);
    end
    [start,fvaltmp] = calcphi2d(im_pair,pupildata.H0,start,zf);
    Zcoefs(indx,:) = start;
    fval(indx) = fvaltmp(end);
  end

  %% If only an image pair is provided, throw out the info on the baseIndex
  if (K == 2)
    fval(params.baseIndex) = [];
    Zcoefs(params.baseIndex,:) = [];
  end
end