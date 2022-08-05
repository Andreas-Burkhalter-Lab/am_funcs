function [u,err,rmg_params] = register_multigrid_auto(u,imM,lambda,rmg_params,tol)
% REGISTER_MULTIGRID_AUTO: non-rigid image registration using a multigrid algorithm
%
% This algorithm performs registration (alignment) of arbitrary-dimensional
% images (typically, 2- or 3-dimensional). It allows "warping," i.e., a
% non-rigid deformation. This function calls REGISTER_MULTIGRID_VCYCLE
% repeatedly until convergence.
%
% Syntax:
%   [u,err,rmg_params] = register_multigrid_auto(uin,imM,lambda,rmg_params,tol)
% where
%   the first 4 inputs and all the outputs are as described in
%     REGISTER_MULTIGRID_VCYCLE;
%   tol (default 0.05) is the amount of fractional decrement we need on
%     each multigrid coarse-grid-correction before deciding that it's not
%     helping anymore
%
% See also: REGISTER_MULTIGRID_VCYCLE.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 5)
    tol = 0.01;
  end

  %% Preserve some information about the input state of rmg_params, so that
  % we can restore it at the end
  keep_n_grids = isfield(rmg_params,'n_grids');
  cg_fval_tol = rmg_params.cg_fval_tol;
  [rmg_params.image_grid.tCumulative] = deal(0);
  
  %% Do the initial v-cycle
  [u,errCur,rmg_params] = register_multigrid_vcycle(u,imM,lambda,rmg_params);
  
  %% Do multiple v-cycles, monitoring the error decrement to figure out how
  % many grids we should use
  err0 = rmg_params.err0;
  errLast = err0;
  err = [errLast errCur];
  n_grids = Inf;
  tCumulative = [rmg_params.image_grid.tCumulative];
  while (n_grids > 1 && 1-errCur/errLast > tol)
    % Do another vcycle, but let's truncate the grids at the point where
    % it stopped making a difference.
    decCGC = [rmg_params.image_grid.err_ratio_decCGC];
    n_grids = find(decCGC < tol,1,'first');
    if isempty(n_grids)
      n_grids = length([rmg_params.image_grid.err_ratio_decCGC])+1;
    end
    if (n_grids < 2)
      warning('register:multigrid','Coarse-grid correction is no longer helping; switching to conjugate gradient');
      n_grids = 1;
      rmg_params.cg_fval_tol = tol/4;
    end
    rmg_params.n_grids = n_grids;
    if ~strcmp(rmg_params.display,'none')
      fprintf('\n\nSetting # of grids to %d\n',n_grids);
    end
    errLast = errCur;
    [rmg_params.image_grid.tCumulative] = deal(0);
    [u,errCur,rmg_params] = register_multigrid_vcycle(u,imM,lambda,rmg_params);
    tCumulative = tCumulative + [rmg_params.image_grid.tCumulative];
    errCur = errCur(end);
    err(end+1) = errCur; %#ok<AGROW>
  end
  
  %% Fix up rmg_params outputs
  if ~keep_n_grids && isfield(rmg_params,'n_grids')
    rmg_params = rmfield(rmg_params,'n_grids');
  end
  rmg_params.err0 = err0;
  rmg_params.cg_fval_tol = cg_fval_tol;
  for i = 1:length(rmg_params.image_grid)
    rmg_params.image_grid(i).tCumulative = tCumulative(i);
  end