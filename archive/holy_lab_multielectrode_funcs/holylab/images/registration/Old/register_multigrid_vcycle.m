function [g,imM,err,rmg_params] = register_multigrid_vcycle(g0,imMoving,lambda,rmg_params)
  down_relaxations = rmg_params.relaxation_iterations(1);
  up_relaxations = rmg_params.relaxation_iterations(2);
  n_grids = length(rmg_params.image_grid);
  if rmg_params.covariant
    warp_options = struct('covariant',true,'sqrt',true);
  else
    warp_options = struct('covariant',false,'sqrt',false);
  end

  %
  % Move from fine to coarse grids
  % 
  imM0 = imMoving;
  g = g0;
  for gridIndex = 1:n_grids-1
    fprintf('\n\ngridIndex = %d\n',gridIndex);
    % Save the current warped image (needed on the way back to finer grids)
    rmg_params.image_grid(gridIndex).imM = imM0;
    % Relax the estimate of g on this grid, and compute the resulting
    % warped image
    if down_relaxations
      psiM_0 = register_im2psi(imM0,rmg_params);
      rmg_params.image_grid(gridIndex).startval = [];  % clear it so initial val is recalculated
      for relaxIndex = 1:down_relaxations
        [g,psiM,sqrtdetJ,err,rmg_params.image_grid(gridIndex)] = ...
          register_multigrid_relax(g,psiM_0,lambda, ...
          rmg_params.covariant,rmg_params.image_grid(gridIndex));
      end
      imM = register_psi2im(psiM,rmg_params);
    else
      imM = register_warp(imM0,g,rmg_params);
    end
    % Save the current status of the deformation
    rmg_params.image_grid(gridIndex).g = g;
    % Restrict the image to the next coarser grid; also initialize the
    % deformation on that grid to be free of warp
    decimate = rmg_params.image_grid(gridIndex+1).decimate;
    imM0 = imrestrict(imM,decimate);
    rmg_params.image_grid(gridIndex+1).sqrtdetJ0 = ...
      imrestrict(sqrtdetJ,decimate);
    fprintf('restrict (final mu = %g)\n',rmg_params.image_grid(gridIndex).mu);
    g = cell(1,rmg_params.n_dims);
  end
  
  %
  % Do the full minimization on the coarsest grid
  %
  fprintf('\n\nFull minimization (gridIndex %d)\n',n_grids);
  [g,imM,rmg_params.image_grid(n_grids)] = ...
      register_multigrid_minimize(g,imM0,lambda,rmg_params.covariant,...
      rmg_params.image_grid(n_grids),rmg_params.max_coarsest_iter);
  
  %
  % Work back up to fine grids, composing the g found on each level with
  % the previous g, and then doing further relaxations
  %
  g_correction = cell(1,rmg_params.n_dims);
  for gridIndex = n_grids-1:-1:1
    fprintf('\n\ngridIndex = %d\n',gridIndex);
    fprintf('Prolong\n');
%     for dimIndex = 1:rmg_params.n_dims
%       g_correction{dimIndex} = register_prolong_g(g{dimIndex}, ...
%         rmg_params.image_grid(gridIndex).sz);
%     end
    g_correction = register_prolong_g(g,rmg_params.image_grid(gridIndex).sz);
    %g_correction = register_g0(size(g_correction{1}));
    g = register_composeg(rmg_params.image_grid(gridIndex).g,...
      g_correction);
    % Test whether this new g is better than the old one (sadly, it's not
    % guaranteed)
    psiM_0 = register_im2psi(rmg_params.image_grid(gridIndex).imM,rmg_params);
    warp_options.sqrtdetJ0 = rmg_params.image_grid(gridIndex).sqrtdetJ0;
    [psig,w,sqrtdetJ] = register_warp(psiM_0,g,warp_options);
    E = register_E(rmg_params.image_grid(gridIndex).psiF,psig,w,sqrtdetJ,...
      rmg_params.image_grid(gridIndex).pixel_spacing,lambda);
    if (E > rmg_params.image_grid(gridIndex).startval)
      % It's worse, go back to the original and use gradient
      % smoothing instead of restriction to calculate the low-frequency
      % corrections
      g = rmg_params.image_grid(gridIndex).g;
      [g,rmg_params.image_grid(gridIndex)] = ...
        register_multigrid_smgrad(g,psiM_0,lambda,warp_options,rmg_params.image_grid(gridIndex));
    else
      rmg_params.image_grid(gridIndex).startval = E; % use the new value as the "best yet"
    end
    if up_relaxations
      for relaxIndex = 1:up_relaxations
        [g,psiM,sqrtdetJ,err,rmg_params.image_grid(gridIndex)] = ...
          register_multigrid_relax(g,psiM_0,lambda,rmg_params.covariant,rmg_params.image_grid(gridIndex));
      end
      imM = register_psi2im(psiM,rmg_params);
    end
  end
