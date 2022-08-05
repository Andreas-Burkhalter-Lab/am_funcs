function [g,psig,mu,old_val] = register_graddescent(psi1,psi2,g,lambda,mu,options)
% REGISTER_GRADDESCENT: non-rigid multidimensional image registration
% This function starts from an initial guess and improves it.  This is
% great when your initial guess is good.  However, you can end up being
% trapped in local minima.  REGISTER_MULTIRES provides a multi-resolution
% approach that at least mitigates this problem, and is probably the best
% place to go if you're starting a problem "from scratch."  This
% function, however, is probably the best bet if you need to improve an
% existing estimate.
%
% Syntax:
%   [g,psig,mu,err] = register_graddescent(psi1,psi2,g0,lambda,mu)
%   [g,psig,mu,err] = register_graddescent(psi1,psi2,g0,lambda,mu,options)
% where
%   psi1 is the square-root of image1 (this does covariant registration);
%   psi2 is the square-root of image2 (the "moving" image);
%   g0 is the initial guess for the deformation (see REGISTER_G0);
%   lambda is the regularization parameter;
%   mu is the initial stepsize parameter for the gradient descent;
%   options is a structure which may have some fields (see source code);
% and
%   g is the optimized deformation;
%   psig is the deformed psi2;
%   mu is the final stepsize parameter;
%   err is the registration error (includes only data terms).
%   
% See also: REGISTER_MULTIRES, REGISTER_G0.
  
% Copyright 2006 by Timothy E. Holy
  
  sz = size(psi1);
  n_dims = ndims(psi1);
  if (nargin < 5)
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
    options.step2_thresh = 1e-6;
  end
  if ~isfield(options,'sqrt')
    options.sqrt = true;
  end
  if isempty(g)
    g = register_g0(psi1);
  end
  % Sanity checks
  cls1 = class(psi1);
  cls2 = class(psi2);
  if (cls1(1) == 'u' || cls2(1) == 'u')
    error('Unsigned types not supported, use regular ints')
  end
  for dimIndex = 1:n_dims
%     if ~isequal(size(g{dimIndex}),sz)
%       error('Size of g does not match image (psi1) size');
%     end
    if (isfield(options,'grad_psi1') && length(options.grad_psi1) >= dimIndex)
      if ~isequal(size(options.grad_psi1{dimIndex}),sz)
        error('Size of grad_psi1 does not match psi1 size');
      end
    else
      % Calculate the gradient of psi1
      options.grad_psi1{dimIndex} = deriv(psi1,dimIndex);
    end
  end
  
  colons = repmat({':'},1,n_dims);
  
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
  %
  % The product of psig & psi1 helps put lambda in terms of natural units
  %
  sqr_img = psig .* psi1;
  sqr_scalar = abs(nanmean(sqr_img(:)));
  clear psi1_laplacian psig_laplacian H_img sqr_img
  %
  % Compute the value of the energy
  %
  old_val = register_gradstep_energy(psi1,psig);
  still_descending = true;
  iter = 0;
  while still_descending
    %
    % Compute the gradient with respect to g
    %
    grad_g = register_gradient(psi1,psig,J,detJ,options);
    %[grad_data,grad_reg] = register_gradient(psi1,psig,J,detJ,options);
    %grad_g = grad_data + (lambda*sqr_scalar)*grad_reg;
    % Reduce the gradient, if necessary (i.e., regularization by grid size)
    sz_grad = size(grad_g);
    if ~isequal(sz_grad(2:end),size(g{1}))
      n_decimations = size(options.g_decimate,1);
      for decimationIndex = 1:n_decimations
        grad_g = imreduce(grad_g,[1 options.g_decimate(decimationIndex,:)]);
      end
    end
    
    is_lower = false;
    while (~is_lower && iter < options.iter_max && still_descending)
      %
      % Update g, and calculate the new energy
      %
      for dimIndex = 1:n_dims
        g_new{dimIndex} = g{dimIndex} - (mu/Hscalar)*permute(grad_g(dimIndex,colons{:}),[2:n_dims+1 1]);
      end
      [psig_new,detJ_new,J_new] = register_warp(psi2,g_new,options);
      new_val = register_gradstep_energy(psi1,psig_new);
      is_lower = (new_val <= old_val);
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
      g = g_new;
      psig = psig_new;
      J = J_new;
      detJ = detJ_new;
      old_val = new_val;
      step2 = (mu/Hscalar)^2 * grad_g.^2;
      step2 = sum(step2,1);
      still_descending = any(step2(:) > options.step2_thresh);
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
  %n_pix = numel(psi1);
  % The data penalty
  nanFlag = isnan(psig);
  tmp_im = (psi1-psig).^2;
  energy = mean(tmp_im(~nanFlag(:)));
  % The regularization penalty
  %tmp_im = log(abs(detJ));
  %energy = energy - lambda*mean(tmp_im(:));  % log|detJ|
  %energy = energy + lambda*mean(tmp_im(:).^2);  %(log|detJ|)^2


function [grad_data,grad_reg] = register_gradient(psi1,psig,J,detJ,options)
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
  % Compute the wronskian
  w = zeros([n_dims sz],'single');
  for dimIndex = 1:n_dims
    w(dimIndex,colons{:}) = psig .* grad_psi1{dimIndex} - ...
      psi1 .* grad_psig{dimIndex};
  end
  w(isnan(w)) = 0;  % Get rid of data terms when over the edge
  % Compute the inverse transposed Jacobian, and while at it multiply the
  % wronskian to get the gradient with respect to the data
  [Jti,grad_data] = register_Jtinvw(J,detJ,w);
  if (nargout > 1)
    % Regularization: laplacian of laplacian of g
    for dimIndex1 = 1:n_dims
      lapg = zeros(size(detJ));
      for dimIndex2 = 1:n_dims
        Jtmp = permute(J(dimIndex1,dimIndex2,colons{:}),[3:n_dims+2 1 2]);
        lapg = lapg + deriv(Jtmp,dimIndex2);
      end
      laplapg = zeros(size(detJ));
      for dimIndex2 = 1:n_dims
        laplapg = laplapg + ...
          deriv(deriv(lapg,dimIndex2),dimIndex2);
      end
      grad_reg(dimIndex1,colons{:}) = 2*laplapg;
    end
  end
  return
  % Old, wrong regularization code
  if (nargout > 1)
    lndetJ = 2*log(abs(detJ));
    % Now do the regularization term
    grad_reg = zeros([n_dims sz],class(grad_data));
    % Compute the divergence of each row
    for dimIndex1 = 1:n_dims
      for dimIndex2 = 1:n_dims
        grad_reg(dimIndex1,colons{:}) = squeeze(grad_reg(dimIndex1,colons{:})) + ...
          deriv(squeeze(Jti(dimIndex1,dimIndex2,colons{:})).*lndetJ,dimIndex2);
      end
    end
  end
  % Old, slow code to calculate Jti follows
  Jt = permute(J,[2 1 3:(n_dims+2)]);
  Jti = nan([n_dims n_dims sz],'single');
  grad_data = nan([n_dims sz],'single');
  for pixelIndex = 1:n_pix
    tmp = inv(Jt(:,:,pixelIndex));
    grad_data(:,pixelIndex) = tmp * w(:,pixelIndex);
    Jti(:,:,pixelIndex) = tmp;
  end
