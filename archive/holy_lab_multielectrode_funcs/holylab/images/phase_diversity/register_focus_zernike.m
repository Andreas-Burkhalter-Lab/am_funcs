function [Zcoefs,fval,obj] = register_focus_zernike(im,pupildata)
% REGISTER_FOCUS_ZERNIKE: calculate tip/tilt and defocus for object images
%
% Syntax:
%   [Zcoefs,fval,obj] = register_focus_zernike(im,pupildata)

   sz = size(im);
   sz(end+1:3) = 1;
   K = sz(3);
   
   % Use the variance as a measure of image sharpness (FIXME?)
   v = zeros(1,K);
   for k = 1:K
     tmp = im(:,:,k);
     v(k) = var(tmp(:));
   end
   % Choose the sharpest image as the "base image"
   [maxv,baseIndex] = max(v);
   
   Zindex = [1 2 4];
   params = struct('Zindex', [1 2 4], 'baseIndex', baseIndex);
   [Zcoefs,fval] = register_zernike(im,pupildata,params);
   
   % Now do another optimization where we don't require a particular "base" image
   config.mode = 'Zn block';
   config.Zindex = [1 2 4];
   phimodel = phi_parametrizations(config,pupildata);
   
%    p0 = zeros(K,3);
%    p0(:,3) = 0.1;
%    
%    % Do a coarse first pass
%    phimodel.tol = 1e-4;
%    fvalOld = Inf;
%    [Zcoefs,fval,obj] = calcphi2d(im,pupildata.H0,p0,phimodel);
%    fval = fval(end);
%    while (fval < fvalOld)
%      fvalOld = fval(end);
%      % Determine which images are well-explained and which are
%      % poorly-explained
%      phi = phimodel.param2phi(Zcoefs);
%      [Hk,sk,imc] = pd_forward_model_2d(phi,pupildata.H0,obj);
%      err = zeros(1,K);
%      for k = 1:K
%        tmp = im(:,:,k) - imc(:,:,k);
%        err(k) = sum(tmp(:).^2);
%      end
%      % Do linear interploation to try to get a good starting guess
%      % We take abs(defocus) because the sign is ambiguous
%      pabs = Zcoefs; pabs(:,3) = abs(Zcoefs(:,3));
%      keepFlag = err < median(err);
%      [Zslope,Zstatic] = linregress(find(keepFlag),pabs(keepFlag,:));
%      Z0 = (1:K)' * Zslope + repmat(Zstatic,K,1);
%      [Zcoefs,fval,obj] = calcphi2d(im,pupildata.H0,Z0,phimodel);
%      fval = fval(end);
%    end
%    % Do a more refined optimization
%    phimodel.tol = 1e-6;

   [Zcoefs,fval,obj] = calcphi2d(im,pupildata.H0,Zcoefs,phimodel);
      