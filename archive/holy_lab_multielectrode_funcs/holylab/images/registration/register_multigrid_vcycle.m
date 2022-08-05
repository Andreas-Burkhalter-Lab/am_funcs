function [u,err,rmg_params] = register_multigrid_vcycle(u,imMoving,lambda,rmg_params)
% REGISTER_MULTIGRID_VCYCLE: non-rigid image registration using a multigrid algorithm
%
% This algorithm performs registration (alignment) of arbitrary-dimensional
% images (typically, 2- or 3-dimensional). It allows "warping," i.e., a
% non-rigid deformation. Note this function does either a v-cycle or
% a w-cycle, depending on the setting of the parameter in rmg_params (see
% below).
%
% Syntax:
%   [u,err,rmg_params] = register_multigrid_vcycle(uin,imM,lambda,rmg_params)
% where
%   uin is the initial guess for the deformation (in the absence of a guess,
%     supply the empty matrix []). The moving image will be evaluated at
%     coordinates x' = g(x) = x + u(x), so u encapsulates the deviation
%     from "no deformation".  Note this will probably not be the same size as
%     the image, as it should be consistent with the "gap_data" field of
%     rmg_params (see REGISTER_MULTIGRID_OPTIONS)
%   imM is the "moving image", the one we want to warp
%   lambda is the penalty on volume changes; increasing lambda will force
%     more "rigid" deformations. An "intermediate" value is to set lambda
%     equal to the square of the expected pixel-by-pixel noise.
%   rmg_params is the output of REGISTER_MULTIGRID_OPTIONS (feel free to
%     make any adjustments in the parameters, of course)
%
% and
%   u is the optimized deformation
%   err is a vector containing the registration penalty on successive
%     v-cycles (the first entry is for the initial unregistered image, the
%     last entry is the final value of the registration penalty)
%   rmg_params, on output, contains a great deal of diagnostic information
%     that may or may not be useful.
% To get the registered image, do this:
%   imMg = register_multigrid_warp(imM,u,rmg_params);
%
% Examples (note that these try to be "realistic" by taking a region from
% a larger deformed image, just as you might expect if you camera is
% looking at a shifting scene):
%
% 1. A rigid deformation (note REGISTER_RIGID is the better choice in such cases)
%   im = double(imread('cameraman.tif'));
%   rng = 50:200;
%   dx = 5; dy = -2; f = 0.7;  % f generates a sub-pixel shift along x-coordinate
%   imM = im(rng,rng);
%   imF = f*im(rng+dx,rng+dy) + (1-f) * im(rng+dx+1,rng+dy);
%   % Show the original images
%   figure; imshowrgb(imF,imM); title('Unregistered images')
%   % Set up the registration
%   mask = true(size(imF)); mask = trim_mask_edges(mask,[10 10]);
%   rmg_params = register_multigrid_options(imF,mask);
%   % (To make things work best, now set rmg_params.g_level_gap higher)
%   lambda = 10*mean(imF(:).^2);  % force a semi-rigid deformation
%   [u,err,rmg_params] = register_multigrid_vcycle([],imM,lambda,rmg_params);
%   % Examine the result
%   imMg = register_multigrid_warp(imM,u,rmg_params);
%   figure; imshowrgb(imF,imMg); title('Registered images')
%   figure; register_visualize_u(u,rmg_params); title('Deformation')
%   figure; plot(0:length(err)-1,err); title('error as a function of iteration')
%
%  The small residual errors in this example are due to the fact that the
%  image isn't defined over-the-edge. Note that in this case, since the
%  image data strongly constrain the result, one can set lambda = 0 and
%  still get a deformation that turns out to be rigid over most of the
%  image area (and indeed the registration result turns out to be of even
%  higher quality).
%
% 2. A non-rigid deformation
%   im = double(imread('cameraman.tif'));
%   % Create a random deformation
%   scale = 10;  % max # of pixels in shift
%   sz = size(im);
%   udef = randn([sz 2]);
%   udef = imfilter_gaussian(udef,[10 10 0]);  % smooth the random deformation
%   udef = udef * (scale / max(abs(udef(:)))); % set the scale
%   imM_full = register_multigrid_warp(im,udef);
%   % Snip out a portion
%   rng = 15:240;
%   imM = imM_full(rng,rng);
%   imF = im(rng,rng);
%   figure; imshowrgb(imF,imM); title('Unregistered images')
%   % Set up the registration
%   mask = true(size(imF)); mask = trim_mask_edges(mask,[1 1]*scale+3);
%   rmg_params = register_multigrid_options(imF,mask);
%   lambda = 1;  % allow a very "floppy" registration
%   [u,err,rmg_params] = register_multigrid_vcycle([],imM,lambda,rmg_params);
%   imMg = register_multigrid_warp(imM,u,rmg_params);
%   figure; imshowrgb(imF,imMg); title('Registered images')
%   figure; imshowrgb(imF.*mask,imMg.*mask); title('Registered images, inside mask')
%   figure; register_visualize_u(u,rmg_params); title('Deformation')
%   figure; plot(0:length(err)-1,err); title('error as a function of iteration')
%
% See also: REGISTER_DEMO, REGISTER_MULTIGRID_OPTIONS, REGISTER_RIGID2NONRIGID, REGISTER_MULTIGRID_AUTO.

% Copyright 2010 by Timothy E. Holy

  %% Initialization
  % Determine the "gaps" (difference in multigrid levels) between the
  % image data, the regularization penalty, and the parametrization of
  % the deformation u. 
  gap_data = rmg_params.gap_data;
  gap_regularize = rmg_params.gap_regularize;
  if (gap_regularize > gap_data)
    error('rmg_params.gap_regularize cannot be larger than rmg_params.gap_data');
  end
  % Determine the number of grids to use. By default use as many as
  % allowed by the data.
  n_grids = length(rmg_params.image_grid) - gap_data;
  if (isfield(rmg_params,'n_grids') && ~isempty(rmg_params.n_grids))
    n_grids = min(n_grids,rmg_params.n_grids);
  end
  while (n_grids > 0 && (isempty(rmg_params.image_grid(n_grids).imFixed) || sum(rmg_params.image_grid(n_grids).mask(:)) == 0))
    n_grids = n_grids - 1;
  end
  rmg_params.level_u_max = n_grids + gap_data;
  
  if strcmp(rmg_params.display,'all')
    fprintf('%d grids will be used. The finest u-grid is of size ',n_grids);
    fprintf('%d ',rmg_params.image_grid(1+gap_data).sz);
    fprintf('\nand the coarsest u-grid is of size ');
    fprintf('%d ',rmg_params.image_grid(n_grids+gap_data).sz);
    fprintf('\n');
%    fprintf('(If you want the coarsest grid to be even smaller,\n consider increasing gap_data.)\n');
  end
  % If necessary, compute the maximum number of conjgrad iterations.
  if strcmp(rmg_params.max_coarsest_iter_mode,'auto')
    numel_u_coarsest = prod(rmg_params.image_grid(n_grids+gap_data).sz);
    rmg_params.max_coarsest_iter = 20*numel_u_coarsest;
  end
  % Set defaults
  if ~isfield(rmg_params,'options')
    rmg_params.options = struct;
  end
  rmg_params.options = default(rmg_params.options,'Mfull',false,'hessmax',10000,'zero_nans',true);
  % Clear old history
  for i = 1:length(rmg_params.image_grid)
    rmg_params.image_grid(i).err_ratio_relax = 1;
    rmg_params.image_grid(i).err_ratio_CGC = 1;
    rmg_params.image_grid(i).err_ratio_comp = 1;
    rmg_params.image_grid(i).muCGC = [];
    rmg_params.image_grid(i).muComp = [];
    rmg_params.image_grid(i).tCumulative = 0;
  end
  rmg_params.err0 = nan;

  %% Run multigrid
  [u,err,rmg_params] = rmg_recurse(u,1,imMoving,lambda,1+gap_data,rmg_params);
end

function [u,fval,rmg_params] = rmg_recurse(u,detJprev,imM,lambda,level_u,rmg_params,fval0)
  %% Handles a single level of multigrid, going first to coarser grids and
  % then improving the current grid
  tStart = tic;
  gap_data = rmg_params.gap_data;
  gap_regularize = rmg_params.gap_regularize;
  n_dims = rmg_params.n_dims;
  level_data = level_u - gap_data;
  level_reg = level_u - gap_regularize;
  imF = rmg_params.image_grid(level_data).imFixed;
  mask = rmg_params.image_grid(level_data).mask;
  % Pre-calculate the identity deformations
  if isfield(rmg_params.image_grid,'u0')
    u0 = rmg_params.image_grid(level_u).u0;
    u0reg = rmg_params.image_grid(level_reg).u0;
  else
    u0 = zeros([rmg_params.image_grid(level_u).sz n_dims],class(imF));
    u0reg = zeros([rmg_params.image_grid(level_reg).sz n_dims],class(imF));
  end
  if isfield(rmg_params.image_grid,'g0')
    g0data = rmg_params.image_grid(level_data).g0;
  else
    g0data = g_cell2array(register_g0(size(imF),class(imF)));
  end
  if isfield(rmg_params,'shift')
    restrictFlag = cat(1,rmg_params.image_grid.restrict);
    cdF = cumsum(restrictFlag,1);
    scale = 2.^cdF(level_data,:);
    colons = repmat({':'},1,rmg_params.n_dims);
    for dimIndex = 1:rmg_params.n_dims
      g0data(colons{:},dimIndex) = g0data(colons{:},dimIndex)+rmg_params.shift(dimIndex)/scale(dimIndex);
    end
  end
    
  upfunc = @(u) register_u_prolong_wrapper(u,rmg_params,{level_data,level_reg},g0data,u0reg);
  downfunc = @(gradu) register_gradu_restrict(gradu,rmg_params,level_u);
    
  if ~isempty(u)
    %% We have an initial guess, do compositional improvement
    % Apply the deformation to the image
    options = rmg_params.options;
    options.Mfull = true;
    [imMg,detJnew,fval0] = register_multigrid_penalty(u,detJprev,imM,imF,mask, ...
      lambda,rmg_params.detJvalue,upfunc,[],options);
    if (level_data == 1 && isnan(rmg_params.err0))
      rmg_params.err0 = fval0;
    end
    % Do multigrid on the deformed image
    options = rmg_params.options;
    fval0Comp = register_valgrad(u0,detJnew.*detJprev,imMg,imF,mask, ...
      lambda,rmg_params.detJvalue,upfunc,[],options);
    if isinf(fval0Comp)
      if strcmp(rmg_params.display,'all')
        fprintf('level_data = %d, skipping composition because of NaN (could mask out more on edges)\n', level_data);
      end
      fval = fval0;
      return
    end
    % Clear out any shift value, so it doesn't get added twice
    if isfield(rmg_params,'shift')
      shift = rmg_params.shift;
      rmg_params = rmfield(rmg_params,'shift');
    end
    [uNew,fval,rmg_params] = rmg_recurse([],detJnew.*detJprev,imMg,lambda,level_u,rmg_params,fval0Comp);
    % Compose the correction with the original deformation
    szu = size(u);
    g0 = g_cell2array(register_g0(szu(1:end-1),class(imF)));
    gNew = g0+uNew;
    uTest = uNew + register_composeg(u,gNew);
    du = uTest-u;
    % "Optimize" the use of the compositional improvement
    % We don't have to optimize this to a high degree, because it should
    % apply in its entirety; if not, there is something wrong
    muComp = 2;
    muCompmin = 0.25;
    fval = fval0+1;
    while (~(fval < fval0) && muComp > muCompmin)
      muComp = muComp/2;
      uTest = u + muComp*du;
      fval = register_valgrad(uTest,detJprev,imM,imF,mask, ...
        lambda,rmg_params.detJvalue,upfunc,[],options);
    end
    if (fval < fval0)
      u = uTest;
    else
      muComp = 0;
      fval = fval0;
    end
    if strcmp(rmg_params.display,'all')
      fprintf('level_data = %d, composition muComp = %g, fractional error decrement = %g\n',level_data,muComp,1-fval/fval0);
    end
    rmg_params.image_grid(level_u).err_ratio_comp = rmg_params.image_grid(level_u).err_ratio_comp * (fval/fval0);
    % Relax the combined deformation
    mu = [];
    fval00 = fval;
    for i = 1:rmg_params.n_relaxations
%       [u,fval,mu] = rmg_improve(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,mu,rmg_params.linemin_fval_tol);
      mu = mu/10;
      [u,fval,tmp,mu] = register_improve_hessian(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,mu,options);
    end
    if (rmg_params.n_relaxations > 0 && strcmp(rmg_params.display,'all'))
      fprintf('level_data = %d, post-composition relaxation: fractional error decrement = %g\n',level_data,1-fval/fval00);
    end
    % Do reporting
    rmg_params.image_grid(level_u).err_ratio_relax = rmg_params.image_grid(level_u).err_ratio_relax*(fval/fval00);
    rmg_params.image_grid(level_u).tCumulative = ...
      rmg_params.image_grid(level_u).tCumulative + toc(tStart);
    return
  else
    %% Initial guess is the identity (empty u), so do standard cycle
    if (level_u >= rmg_params.level_u_max)
      % We're on the coarsest grid, so iterate to convergence
      fval = 0;
      fval0 = 1;
      if isempty(u)
        u = u0;
      end
      first = true;
      while ((fval0-fval) > rmg_params.coarsest_fval_tol*(fval0+fval))
        [u,fval,fval0] = register_improve_hessian(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,[],rmg_params.options);
        if first
          fval00 = fval0;
          first = false;
        end
      end
      if strcmp(rmg_params.display,'all')
        fprintf('level_data = %d, relaxation: fractional error decrement = %g\n',level_data,1-fval/fval00);
      end
      rmg_params.image_grid(level_u).err_ratio_relax = rmg_params.image_grid(level_u).err_ratio_relax*(fval/fval00);
      rmg_params.image_grid(level_u).tCumulative = ...
        rmg_params.image_grid(level_u).tCumulative + toc(tStart);
      return
      % We're on the lowest level, use conjugate gradient
%       mu = rmgv_setmu(rmg_params,level_u);
%       cgfunc = @(u) register_valgrad(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc);
%       cgopts = struct('iter_max',rmg_params.max_coarsest_iter,...
% 		      'mu',mu,'tol',rmg_params.cg_fval_tol);
%       [u,fval] = conjgrad(cgfunc,u0,cgopts);
%       [u,fval,exitflag,output] = fminunc(cgfunc,u,optimset('GradObj','on'));
    else
      %% Do a coarse-grid correction
      % Step 1: coarsen the moving image
      restrictFlagIm = rmg_params.image_grid(level_data+1).restrict;
      imMH = array_restrict(imM,restrictFlagIm);
      if (lambda > 0)
        % Restrict the prior determinant penalty
        % Since restriction effectively assumes that beyond-the-edge
        % values are 0, transform detJprev so that 0 corresponds to the
        % value we want (using the same scaling as the penalty)
        dJ = log(detJprev/rmg_params.detJvalue);
        detJH = rmg_detJ_restrict(dJ,rmg_params,level_reg);
        detJH = exp(detJH)*rmg_params.detJvalue;
      else
        detJH = rmg_params.detJvalue;
      end
      % Step 2: solve on coarser grid.
      [uH,err,rmg_params] = rmg_recurse([],detJH,imMH,lambda,level_u+1,rmg_params);
      % Step 3: bring correction up to finer grid
      uh = register_u_prolong(uH,rmg_params,level_u);
      % Step 4: "optimize" the use of the coarse-grid correction
      % We do this quite crudely, because if muCGC is not near 1,
      % coarse-grid isn't really helping like one might hope it would
%       [muCGC,fval,u,fval0] = ...
%         rmg_linesearch(u0,uh,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,1,[],rmg_params.linemin_fval_tol);
      if (nargin < 7)
        fval0 = register_valgrad(u0,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,rmg_params.options);
      end
      if (level_data == 1 && isnan(rmg_params.err0))
        rmg_params.err0 = fval0;
      end
      muCGC = 2;
      muCGCmin = 0.25;
      fval = fval0+1;
      while (~(fval < fval0) && muCGC > muCGCmin)
        muCGC = muCGC/2;
        u = muCGC*uh;
        fval = register_valgrad(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,rmg_params.options);
      end
      if (fval >= fval0)
        u = u0; % give up on coarse-grid correction
        fval = fval0;
        muCGC = 0;
      end
      if ~strcmp(rmg_params.display,'none')
        fprintf('level_data = %d, muCGC = %g, fractional error decrement %g\n',level_data,muCGC,1-fval/fval0);
      end
      rmg_params.image_grid(level_u).muCGC(end+1) = muCGC; % report back how it worked
      rmg_params.image_grid(level_u).err_ratio_CGC = rmg_params.image_grid(level_u).err_ratio_CGC * (fval/fval0);
      % Step 5: perform fine-grid polishing
      fval0 = fval;
      mu = [];
      for i = 1:rmg_params.n_relaxations
        mu = mu/10;
%         [u,fval,mu] = rmg_improve(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,mu,rmg_params.linemin_fval_tol);
        [u,fval,tmp,mu] = register_improve_hessian(u,detJprev,imM,imF,mask,lambda,rmg_params.detJvalue,upfunc,downfunc,mu,rmg_params.options);
      end
      if (rmg_params.n_relaxations > 0 && strcmp(rmg_params.display,'all'))
        fprintf('level_data = %d, relaxation: fractional error decrement = %g\n',level_data,1-fval/fval0);
      end
      rmg_params.image_grid(level_u).mu = mu;
      rmg_params.image_grid(level_u).err_ratio_relax = rmg_params.image_grid(level_u).err_ratio_relax * (fval/fval0);
      rmg_params.image_grid(level_u).tCumulative = ...
        rmg_params.image_grid(level_u).tCumulative + toc(tStart);
      % Step 6: decide whether to do a wcycle
      if (rmg_params.wcycle && level_data > 1)
        rmg_params.wcycle = false;
        [u,err,rmg_params] = rmg_recurse(u,detJprev,imM,lambda,level_u,rmg_params);
        rmg_params.wcycle = true;
      end
    end
  end
end
	

function varargout = register_valgrad(varargin)
  tmp = cell(1,nargout+2);
  varargout = cell(1,nargout);
  [tmp{:}] = register_multigrid_penalty(varargin{:});
  for i = 1:nargout
    varargout{i} = tmp{i+2};
  end
end

% function [imMg,detJnew,v,vgrad] = register_penalty(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,Mfull)
%   if (nargin < 10)
%     % By default, do not evaluate the complete imMg, just evaluate
%     % the pixels we need
%     Mfull = false;
%   end
%   detJnew = 1;
%   tmp = cell(1,nargout-1);
%   if (lambda > 0)
%     [ghdata,uhreg] = upfunc(u);
%   else
%     ghdata = upfunc(u);
%   end
%   [tmp{:}] = register_mismatch_noncov(ghdata,imM,imF,mask,Mfull);
%   imMg = tmp{1};
%   p1 = tmp{2};
%   if (nargout > 3)
%     grad1 = tmp{3};
%   end
%   if (lambda > 0)
%     if isscalar(detJprev)
%       [tmp{:}] = register_logdetpenalty(uhreg,detJvalue/detJprev);
%     else
%       [tmp{:}] = register_logdetpenalty(uhreg,detJvalue,detJprev);
%     end
%     detJnew = tmp{1};
%     p2 = tmp{2};
%     if (nargout > 3)
%       grad2 = tmp{3};
%     end
%   end
%   v = p1;
%   if (lambda > 0)
%     v = v + lambda*p2;
%   end
%   if (nargout > 3)
%     % Bring the grad back down to the size of u
%     grad1 = downfunc(grad1);
%     vgrad = grad1;
%     if (lambda > 0)
%       grad2 = downfunc(grad2);
%       vgrad = vgrad + lambda*grad2;
%     end
%   end
%   if isnan(v)
%     %warning('Mismatch value is NaN; image may be shifted beyond field of view');
%     v = Inf;
%     if (nargout > 3)
%       vgrad = zeros(size(vgrad));
%     end
%   end
%   if (nargout > 3)
%     if any(isnan(vgrad(:)))
%       error('gradient should not have NaNs...');
%     end
%   end
% end

% function [u,fval,mu,fval0] = rmg_improve(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,mu,tolfun)
%   [fval0,fgrad] = register_valgrad(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc);
%   [mu,fval,u] = rmg_linesearch(u,-fgrad,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,mu,fval0,tolfun);
% end

% function [mu,fval,unew,fval0] = rmg_linesearch(u,du,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,mu,fval0,tolfun)
%   linfunc = @(mu) register_valgrad(u+mu*du,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc);
%   if isempty(fval0)
%     fval0 = linfunc(0);
%   end
%   [mutriple,fvaltriple] = linbrackgd(linfunc,fval0,-sum(du(:).^2),mu);
%   if (length(mutriple) > 1)
%     [mu,fval] = brentmin(linfunc,mutriple,fvaltriple,struct('TolFun',tolfun));
%   else
%     mu = 0;
%     fval = fval0;
%   end
%   unew = u+mu*du;
% end

function [unew,v,v0,mu] = register_improve_hessian(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,mu,options)
  if (nargin < 10 || isempty(mu))
    mu = 1e-2;
  end
  mumax = 1e6;  % keep it from iterating forever when no improvement is possible
  [v0,vgrad,hess] = register_valgrad(u,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,options);
  if ~isfinite(v0)
    error('register:hessian','The original value must be finite, aborting. Note sum(mask(:)) = %g\n',sum(mask(:)));
  end
  n = size(hess,1);
  hessdiag = spdiags(hess,0);
  meandiag = mean(hessdiag);
  hessdiag = spdiags(meandiag(ones(n,1)),0,n,n);
  v = v0+1;
  while ~(v < v0) && mu < mumax
    mu = mu*10;
    if (numel(vgrad) < options.hessmax)
      % Do Newton-type step
      hessnew = hess + mu*hessdiag;
      du = (-sqrt(1+mu))*(hessnew \ vgrad(:)); % scaling lets it switch to grad descent without making a tiny step
      du = reshape(du,size(vgrad));
    else
      % Hessian is too big, the inversion will take too long, just do it by
      % gradient descent
      du = vgrad/(-meandiag*sqrt(1+mu));
    end
    if ~isempty(u)
      unew = u+du;
    else
      unew = du;
    end
    v = register_valgrad(unew,detJprev,imM,imF,mask,lambda,detJvalue,upfunc,downfunc,options);
  end
end


function dJ = rmg_detJ_restrict(dJ,rmg_params,level_reg)
  % detJ is computed on half-grid points, so we need to make
  % appropriate adjustments to restriction: for dimensions in which u has
  % an odd number of points, we first embed in a larger array, and
  % after restriction the edges are trimmed
  if isscalar(dJ)
    return
  end
  lvlIndex = level_reg+1;
  n_dims = rmg_params.n_dims;
  dJsz = size(dJ);
  restrictFlagU = rmg_params.image_grid(lvlIndex).restrict;
  isOdd = mod(dJsz-1,2) & restrictFlagU;  % is _u_ odd
  if any(isOdd)
    x = cell(1,n_dims);
    dJtmp = dJ;
    dJ = zeros(dJsz+2*isOdd);
    for dimIndex2 = 1:n_dims
      x{dimIndex2} = (1:dJsz(dimIndex2)) + isOdd(dimIndex2);
    end
    dJ(x{:}) = dJtmp;
  end
  dJ = array_restrict(dJ,restrictFlagU);
  if any(isOdd)
    for dimIndex2 = 1:n_dims
      if isOdd(dimIndex2)
        x{dimIndex2} = 2:size(dJ,dimIndex2)-1;
      else
        x{dimIndex2} = ':';
      end
    end
    dJ = dJ(x{:});
  end
end

function g = g_cell2array(gc)
  n_dims = length(gc);
  g = cat(n_dims+1,gc{:});
end

function [ghdata,uhreg] = register_u_prolong_wrapper(u,rmg_params,lvl,g0data,u0reg)
  if ~isempty(u)
    if (nargout > 1)
      [uhdata,uhreg] = register_u_prolong(u,rmg_params,lvl);
    else
      uhdata = register_u_prolong(u,rmg_params,lvl{1});
    end
    ghdata = uhdata + g0data;
  else
    ghdata = g0data;
    if (nargout > 1)
      uhreg = u0reg;
    end
  end
end

% function mu = rmgv_setmu(rmg_params,level_u)
%   % Load the default linesearch step size
%   mu = rmg_params.image_grid(level_u).mu;
%   if (isempty(mu) && level_u < length(rmg_params.image_grid))
%     mu = rmg_params.image_grid(level_u+1).mu/10; % initialize from coarser level
%   end
%   if (isempty(mu) || mu == 0)
%     mu = 1;   % the option of last resort!
%   end
% end
%   