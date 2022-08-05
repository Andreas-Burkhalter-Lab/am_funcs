function [ynew_best,neighborIndex_best,n_best,scale_best] = vmams_step1(x,y,options)
% VMAMS_STEP1: move a landmark by variable metric adaptive mean shift
%
% This function varies the scaling of the coordinates until the neighbors
% of the landmark, as judged by MSAMS, are "isotropic" (meaning the
% variance is unity in each coordinate).
%
% Syntax:
%   ynew = vmams_step1(x,y)
%   [ynew,neighborIndex,n] = vmams_step1(x,y)
%
% See MSAMS_STEP1 for an explanation of these inputs and outputs.
%
%
%   [ynew,neighborIndex,n,scale] = vmams_step1(x,y,options)
%
% options is a structure. It gets passed on to msams_step1, so all of
% those options are applicable here. The additional fields are:
%   scale (default 1): the initial scaling applied to the coordinates;
%   anisotropy_thresh (default 1e-8): a convergence tolerance, the allowed
%     discrepancy from unity variance in each coordinate
%   max_iter_vms (default 100): how many times to try for unity variance
%     before giving up.
% scale is the optimal scaling applied to each coordinate.  Note that one
%   can calculate the equivalent of the "msd" output of MSAMS_STEP1 as
%   1./scale.^2.
%
% See also: MSAMS, MSAMS_STEP1.

% Copyright 2006 by Timothy E. Holy
  
  debug_fig = false;
  debug_text = false;
  if debug_fig
    figure
  end
  d = size(x,1);
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'max_iter_vms')
    max_iter_vms = 100;
  else
    max_iter_vms = options.max_iter_vms;
  end
  if ~isfield(options,'anisotropy_thresh')
    %options.anisotropy_thresh = 1e-8;
    options.anisotropy_thresh = 0.01;
  end
  if isfield(options,'scale')
    scale = options.scale(:);
  else
    scale = ones(d,1);
  end
  if ~isfield(options,'robust')
    options.robust = false;
  end
  if ~isfield(options,'min_n_to_change_scale')
    options.min_n_to_change_scale = 0;
  end
%   if ~isfield(options,'min_to_check')
%     options.min_to_check = 30;
%   end
  
  isdone = false;
  iter = 0;
  scale_history = zeros(d,0);  % will use to detect cycling
  msd_dev_best = inf(d,1);
%  msdPrev = [];
  % Try MSAMS_STEP1 until the contributing points have isotropic variance
  while ~isdone
    options.scale = scale;
    if debug_text
      fprintf('  %g %g',scale);
    end
    scale_history(:,end+1) = scale;
    [ynew,neighborIndex,ncur,msd] = msams_step1(x,y,options);
%     if (length(indexContrib) < options.min_n_to_change_scale)
%       %scale_best = scale;
%       return
%     end
    n_scale = max(ncur,options.min_n_to_change_scale);
%     if (~isempty(msdPrev) && any(abs(log(msd) - log(msdPrev)) > 1))
%       'big change'
%     end
%     msdPrev = msd;
    if debug_fig
      plot(x(1,:),x(2,:),'.'); hold on; plot(x(1,neighborIndex(1:ncur)),x(2,neighborIndex(1:ncur)),'g.'); plot(ynew(1),ynew(2),'rx'); hold off; title('vm')
      pause
    end
    if options.robust
      % Recompute a robust analog of the msd
      dy = x(:,neighborIndex(1:n_scale)) - repmat(y,[1 n_scale]);
      dy = dy .* repmat(scale,[1 n_scale]);
      msd = mean(abs(dy),2);
    end
    if debug_text
      fprintf(', %d, %g %g\n',ncur,msd);
    end
    msd_dev = abs(msd-1); % needed for convergence testing
    msd_dev(msd == 0) = 0; % For coordinates with no variance
    if (sum(msd_dev) < sum(msd_dev_best))
      % Save the best yet, in case convergence is not monotonic
      msd_dev_best = msd_dev;
      scale_best = scale;
      ynew_best = ynew;
      neighborIndex_best = neighborIndex;
      n_best = ncur;
    end
    % Now change the scaling to make contributing points more isotropic
    msd(msd == 0) = 1;  % avoid divide by zero
    if options.robust
      scale = scale ./ msd;
    else
      scale = scale ./ sqrt(msd);
    end
    iter = iter+1;
    %isdone = all(msd_dev < 1/n_scale) || ...
    isdone = all(msd_dev < options.anisotropy_thresh) || ...
      any(all(scale_history == repmat(scale,1,size(scale_history,2)),1)) || ...
      iter >= max_iter_vms;
  end
  if debug_text
    fprintf('done\n');
  end
  if (iter >= max_iter_vms)
    warning('VMAMS:noConvergence','VMAMS failed to converge to isotropic variance');
  end
  if debug_fig
    title('Done!')
    pause
    close(gcf)
  end
  