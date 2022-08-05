function [g,psig,mu,old_val] = register_smgrad(psi1,psi2,g,smooth_size, ...
					       mu,options)
% REGISTER_SMGRAD: registration by smoothed-gradient descent
% Given two input (multidimensional) images, calculate the deformation
% that optimizes their alignment (registration).  The problem is
% under-determined without constraints on the deformation, which here is
% set to be smooth by virtue of the descent algorithm, which follows the
% smoothed gradient of the mismatch between the images.
%
% Syntax:
%   [g,psig,mu,err] =
%      register_smgrad(psi1,psi2,g,sigma,mu,options)
% Here,
%   psi1 = sqrt(im1), psi2 = sqrt(im2);
%   g is the initial guess for the deformation (leave blank if you want
%     to start with the identity), see REGISTER_G0;
%   sigma is a vector giving the smoothing length scale (in pixels) in
%     each dimension of the image;
%   mu is a guess for distance to follow the gradient (usually set to 1);
%   options is an optional structure that may have the following fields:
%     scale: the rate of relaxation in each coordinate (default [1 1 1]
%       or however many dimensions your images have...);
%     covariant: if true, the mismatch penalty is covariant under
%       deformations (default true);
%     sqrt: if true, assumes psi1 and psi2 are the square roots of
%       images, otherwise assumes that the first two inputs are the
%       images themselves (default true);
%     iter_max: the maximum number of iterations (default 50);
%     step2_thresh: the threshold on the square step size for determining
%       that no further progress is being made (default 1e-4);
%     mu_increase_factor: the factor by which to multiply mu after a
%       successful step (default 2);
%     mu_decrease_factor: the factor by which to divide mu after a
%       failed step (default 5).
%
% There is a variation on this calling syntax, used when you want to use
% the FFT to peform the smoothing:
%   [g,psig,mu,err] =
%      register_smgrad(psi1,psi2,g,fftinfo,mu,options)
% where the fftinfo input comes from REGISTER_FFTPREP.
%
% On output,
%   g is the optimized deformation;
%   psig is the warped image;
%   mu is the last value of mu used (may not be useful for future steps);
%   err is the mismatch between the registered images.
%
% See also: REGISTER_RUN, REGISTER_FFTPREP, REGISTER_NONRIGID.

% Copyright 2006 by Timothy E. Holy

  sz = size(psi1);
  n_dims = ndims(psi1);
  if (nargin < 6)
    options = struct;
  end
  if ~isfield(options,'mu_decrease_factor')
    options.mu_decrease_factor = 5;
  end
  if ~isfield(options,'mu_increase_factor')
    options.mu_increase_factor = 2;
  end
  if ~isfield(options,'iter_max')
    options.iter_max = 50;
  end
  if ~isfield(options,'step2_thresh')
    options.step2_thresh = 1e-4;
  end
  if ~isfield(options,'sqrt')
    options.sqrt = true;
  end
  if ~isfield(options,'covariant')
    options.covariant = true;
  end
  if ~isfield(options,'use_fft')
    options.use_fft = isstruct(smooth_size);
  end
  if ~isfield(options,'scale')
    options.scale = ones(1,n_dims);
  end
  
  if isempty(g)
    g = register_g0(size(psi1));
  end
  % Sanity checks
  cls1 = class(psi1);
  cls2 = class(psi2);
  if (cls1(1) == 'u' || cls2(1) == 'u')
    error('Unsigned types not supported, use regular ints')
  end
  if isempty(g)
    g = register_g0(sz);
  end
  for dimIndex = 1:n_dims
    if ~isequal(size(g{dimIndex}),sz)
      error('Size of g does not match image (psi1) size');
    end
    if (isfield(options,'grad_psi1') && length(options.grad_psi1) >= dimIndex)
      if ~isequal(size(options.grad_psi1{dimIndex}),sz)
        error('Size of grad_psi1 does not match psi1 size');
      end
    else
      % Calculate the gradient of psi1
      options.grad_psi1{dimIndex} = deriv(psi1,dimIndex);
    end
  end
  %options.output_size = sz;
  
  colons = repmat({':'},1,n_dims);
  
  % Prepare the smoothing schedule to accomodate the fact that, for the
  % gradient, the first dimension is the gradient coordinate and should
  % not be smoothed over
  %n_smooth = size(smooth_schedule,1);
  %smooth_schedule_pad = [zeros(n_smooth,1) smooth_schedule];
  
  %
  % Compute the deformed image
  %
  [psig,detJ,J] = register_warp(psi2,g,options);
  %
  % Estimate the "Hessian" as a single scalar
  % This will help us set the stepsize in natural units
  %
  psi1_laplacian = zeros(size(psi1),class(psi1));
  psig_laplacian = zeros(size(psig),class(psig));
  for dimIndex = 1:n_dims
    psi1_laplacian = psi1_laplacian + ...
      deriv(deriv(psi1,dimIndex),dimIndex);
    psig_laplacian = psig_laplacian + ...
      deriv(deriv(psig,dimIndex),dimIndex);
  end
  H_img = psig .* psi1_laplacian + psi1 .* psig_laplacian;
  Hscalar = abs(nanmean(H_img(:)));
  if (Hscalar == 0)
    warning('It may be that no useful data are available');
  end
  clear psi1_laplacian psig_laplacian H_img
  %
  % Compute the value of the energy
  %
  old_val = register_gradstep_energy(psi1,psig);
  still_descending = true;
  iter = 0;
  while still_descending
    %
    % Compute the gradient of the data terms with respect to g
    %
    grad_g = register_gradient(psi1,psig,J,detJ,options);
    grad_g = repmat(options.scale(:),[1 size(psig)]) .* grad_g;
    %
    % Smooth grad_g
    %
    %grad_g = imfilter_recursive(grad_g,smooth_schedule_pad);
    grad_g_old = grad_g;
    if options.use_fft
      for dimIndex = 1:n_dims
        ggtmp = squeeze(grad_g(dimIndex,smooth_size.padcoords{:}));
        %ggtmp = imdilatenan(ggtmp);
        ggtmp(isnan(ggtmp)) = 0;
        ggfft = fftn(ggtmp);
        ggtmp = ifftn(ggfft .* smooth_size.hfft);
        grad_g(dimIndex,colons{:}) = ggtmp(smooth_size.snipcoords{:});
      end
    else
      %grad_g = imfilter_gaussian(grad_g,[0 smooth_size],'iir');
      grad_g = imfilter_gaussian(grad_g,[0 smooth_size]);
      % Replace all the NaNs with finite values (these arise because
      % registration can introduce NaNs from out-of-field pixels)
      for dimIndex = 1:n_dims
        grad_g(dimIndex,colons{:}) = imdilatenan(grad_g(dimIndex,colons{:}));
      end
    end
    nanflag = isnan(grad_g);
    nnan = sum(nanflag(:));
    is_lower = false;
    while (~is_lower && iter < options.iter_max && still_descending)
      %
      % Update g by moving against the gradient, and calculate the new energy
      %
      delta_g = (-mu/Hscalar) * grad_g;
      for dimIndex = 1:n_dims
        g_new{dimIndex} = g{dimIndex} + permute(delta_g(dimIndex,colons{:}),[2:n_dims+1 1]);
      end
      [psig_new,detJ_new,J_new] = register_warp(psi2,g_new,options);
      new_val = register_gradstep_energy(psi1,psig_new);
      is_lower = (new_val <= old_val);
      [mu old_val new_val]
      if ~is_lower
        mu = mu/options.mu_decrease_factor;
        step2 = (mu/Hscalar)^2 * grad_g.^2;
        step2 = sum(step2,1);
        still_descending = any(step2(:) > options.step2_thresh);
      end
      iter = iter+1;
    end
    if is_lower
      %
      % We were successful, keep the step & determine if converged
      %
      for dimIndex = 1:length(g)
        dgtmp = g_new{dimIndex} - g{dimIndex};
        fprintf('%g ',nanmean(dgtmp(:)));
      end
      fprintf('\n');
      g = g_new;
      psig = psig_new;
      J = J_new;
      detJ = detJ_new;
      old_val = new_val;
      %step2 = (mu/Hscalar)^2 * grad_g.^2;
      %step2 = sum(step2,1);
      %still_descending = any(step2(:) > options.step2_thresh);
      mu = mu * options.mu_increase_factor;
    else
      %
      % We failed to update successfully, we must have run out of
      % iterations
      %
      if (iter >= options.iter_max)
        warning(['Maximum number of iterations exceeded without convergence,' ...
          ' aborting']);
        still_descending = false;
      end
    end
  end  % while still_descending
  
  
function energy = register_gradstep_energy(psi1,psig)
  nanFlag = isnan(psig) | isnan(psi1);
  tmp_im = (psi1-psig).^2;
  if (sum(~nanFlag(:)) > 0)
      energy = mean(tmp_im(~nanFlag(:)));
  else
      energy = inf;
  end
  
    
    
function grad_data = register_gradient(psi1,psig,J,detJ,options)
  if (nargin < 4)
    options = struct;
  end
  n_dims = ndims(psi1);
  sz = size(psi1);
  if ~isequal(size(psig),sz)
    error('psi1 and psig must be of the same size');
  end
  %n_pix = numel(psi1);
  colons = repmat({':'},1,n_dims);
  % Compute gradients
  grad_psi1 = options.grad_psi1;
  grad_psig = cell(1,n_dims);
  for dimIndex = 1:n_dims
    grad_psig{dimIndex} = deriv(psig,dimIndex);
  end
  % Compute the "wronskian"
  w = zeros([n_dims sz],'single');
  if options.covariant
    for dimIndex = 1:n_dims
      w(dimIndex,colons{:}) = psig .* grad_psi1{dimIndex} - ...
	  psi1 .* grad_psig{dimIndex};
    end
  else
    for dimIndex = 1:n_dims
      w(dimIndex,colons{:}) = 2*(psig - psi1) .* grad_psig{dimIndex};
    end
  end
  %w(isnan(w)) = 0;  % Get rid of data terms when over the edge
		    % Compute the inverse transposed Jacobian, and
		    % while at it multiply the
		    % wronskian to get the gradient with respect to the data
  grad_data = register_Jtinvw(J,detJ,w);
