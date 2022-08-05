function [dxOut,dx2Out,n,n_z,sq_summed_disp,summed_sq_disp,options] = bn_preordered(dx,options)
% BN_PREORDERED: the balanced neighborhood criterion
%
% This is the "lowest-level" function, which implements the balanced
% neighborhood criterion. It is assumed that you have a set of points that
% are candidates for inclusion in the neighborhood, and you know the order
% in which you want them considered.  Note that there are other
% higher-level functions (e.g., BN_BY_DISTANCE) that will order the points
% according to some criterion and call this function.
%
% Syntax:
%   [dxOut,dx2Out,n] = bn_preordered(dx)
%   [dxOut,dx2Out,n] = bn_preordered(dx,options)
%   [dxOut,dx2Out,n,n_z,sq_summed_disp,summed_sq_disp,options] = bn_preordered(dx,...)
% where
%   dx is a d-by-N matrix, each column containing the displacement of a
%     data point from the current probe point.  These points are to be
%     supplied in the order in which you want them considered for inclusion
%     in the neighborhood.
% and
%   dxOut is d-by-1, the mean displacement over the balanced neighborhood
%   dx2Out is d-by-1, the mean square displacement over the balanced neighborhood
%   n is the number of points in the balanced neighborhood
%   n_z is the number of points that were needed to exceed the
%     balancing criterion by a factor of options.z
%   sq_summed_disp and summed_sq_disp are the square summed displacement and
%     summed square displacement, respectively, used in the balancing
%     criterion
%
% The options structure may have the following fields:
%   z (default 3): the factor by which the mean displacement must
%     exceed the uncertainty in the displacement to terminate the growth of
%     the neighborhood. This is the analog of a "p-value" (really, a
%     z-score) expressing your threshold for confidence that the mean of
%     the neighborhood has been displaced by more than the noise.
%   backtrack ('fixed','triangle','dprime','outlier', or 'none'): if
%     anything other than 'none', the last few points (typically, O(z^2) of
%     them) will be discarded from the neighborhood. This prevents
%     contamination from outliers, and is recommended.
%   n_min (default 0): the minimum number of points that can be considered
%     to form a neighborhood. 0 would indicate that there are no points in
%     the neighborhood of the probe point (this outcome is possible only if
%     backtrack is true).
%   plot_criterion (default false): if true, plots the data used to
%     evaluate the balancing criterion.  Lines terminated with an 'x'
%     indicate that the criterion was satisfied before running out of
%     points.
% In principle, you should not need to adjust these settings (with the
% possible exception of 'z' to enforce a different p-value).
%
% Example:
%    dxbase = [zeros(1,10) zeros(1,20)+3];
%    dx = dxbase + randn(size(dxbase));
%    [dxOut,dx2Out,n] = bn_preordered(dx);
% It is likely that you will obtain a value of n near 10.
%
% See also: BN_BY_DISTANCE.

% Copyright 2009 by Timothy E. Holy

  [d,N] = size(dx);
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'z',3,'backtrack','fixed','n_min',0,'plot_criterion',false);
  options.n_min = min(options.n_min,N);
  isTriangle = strcmp(options.backtrack,'triangle');
  
  zs = zeros(d,1);
  dxsum = zeros(d,1);
  dx2sum = zeros(d,1);
  z2 = options.z^2;
  if isTriangle
    L = ceil(2*z2);
    nL = min(N,L);
    w = 1 - (1:nL)/L;
    wrep = repmat(w,d,1);
    dxWsum = sum(wrep.*dx(:,1:nL),2);
    dx2Wsum = sum(sum(wrep.*dx(:,1:nL).^2,2));
    dxSsum = sum(dx(:,1:nL),2);
    dx2Ssum = sum(sum(dx(:,1:nL).^2,2));
    sq_summed_disp = sum(dxWsum.^2);
    if (sq_summed_disp > z2 * dx2Wsum)
      dxOut = zs;
      dx2Out = zs;
      n = 0;
      n_z = L;
      summed_sq_disp = dx2Wsum;
      return
    end
  end
  sq_summed_disp = zeros(1,N);
  summed_sq_disp = zeros(1,N);
  for n = 1:N
    thisdx = dx(:,n);
    dxsum = dxsum + thisdx;
    dx2sum = dx2sum + thisdx.^2;
    if isTriangle
      dxWsum = dxWsum + dxSsum/L;
      dx2Wsum = dx2Wsum + dx2Ssum/L;
      dxSsum = dxSsum - thisdx;
      dx2Ssum = dx2Ssum - sum(thisdx.^2);
      if (n+L <= N)
        nextdx = dx(:,n+L);
        dxSsum = dxSsum + nextdx;
        dx2Ssum = dx2Ssum + sum(nextdx.^2);
      end
      dxsum2 = sum(dxWsum.^2);
      dx2sum_tot = dx2Wsum;
    else
      dxsum2 = sum(dxsum.^2);
      dx2sum_tot = sum(dx2sum);
    end      
    sq_summed_disp(n) = dxsum2;
    summed_sq_disp(n) = dx2sum_tot;
    satisfied = dxsum2 > z2*dx2sum_tot; % || sum(thisdx.^2) > z2*dx2sum_tot/n;
    if satisfied
      break
    end
  end
  n_z = n;
  sq_summed_disp = sq_summed_disp(1:n_z);
  summed_sq_disp = summed_sq_disp(1:n_z);
  % We will exclude some of the last points added, because almost by
  % definition they are outliers (after all, they were able to bias the
  % mean beyond statistical significance).
  if satisfied
    switch options.backtrack
      case 'triangle'
        n_z = n+L;
        n = ceil(n/2);
        n = max(n,options.n_min);
        dxsum = sum(dx(:,1:n),2);
        dx2sum = sum(dx(:,1:n).^2,2);
      case 'fixed'
        n = ceil((n_z-z2)/2);
        n = max(n,options.n_min);
        dxsum = sum(dx(:,1:n),2);
        dx2sum = sum(dx(:,1:n).^2,2);
      case 'dprime'
        % "Peel off" points in reverse order, defining a "core" and a
        % "shell."  Maximize the discriminability between these.
        dxsum_shell = zeros(d,1);
        dx2sum_shell = zeros(d,1);
        dprime2 = zeros(1,n_z+1); % 1st entry is n_core=0
        for n_core = n-1:-1:options.n_min
          thisdx = dx(:,n_core+1);
          % Move a point from the core to the shell
          dxsum = dxsum - thisdx;
          dx2sum = dx2sum - thisdx.^2;
          dxsum_shell = dxsum_shell + thisdx;
          dx2sum_shell = dx2sum_shell + thisdx.^2;
          n_shell = n_z-n_core;
          % Determine how well the shell can be discrimiated from the core
          mu_shell = dxsum_shell/n_shell;
          var_shell = dx2sum_shell/n_shell - mu_shell.^2;
          n_core_reg = max(n_core,1);
          mu_core = dxsum/n_core_reg;
          var_core = dx2sum/n_core_reg - mu_core.^2;
          if (n_shell == 1)
            var_shell = var_core;
          end
          dprime2(n_core+1) = sum((mu_core-mu_shell).^2) / (sum(var_core)+sum(var_shell)+eps);
        end
        dprime2_thresh = z2 ./ min([1 1:n_z], [n_z n_z:-1:1]);
        dprime2(dprime2 < dprime2_thresh) = 0;
        [dpmax,n] = max(dprime2);
        n = n-1;  % to compensate for the n_core=0 term
        if (dpmax == 0)
          n = n_z;
        end
        %n = max(n,options.n_min);
        %       if (dpmax < z2)
        %         n = n_z;
        %       end
        dxsum = sum(dx(:,1:n),2);
        dx2sum = sum(dx(:,1:n).^2,2);
      case 'outlier'
        % We fit to an outlier model in which the last k samples contribute to
        % the mean and variance as if they were outliers
        discrep = (n_z-1:-1:0) .* (summed_sq_disp(n_z)-summed_sq_disp) ...
          - (sq_summed_disp(n_z)-sq_summed_disp);
        offset = floor(z2);
        [mindiscrep,n] = min(abs(discrep(1:n_z-offset)));
        if isempty(n)
          n = 0;
        else
          n = n-1;
        end
        n = max(n,options.n_min);
        % Correct dxsum and dx2sum
        dxkill = dx(:,n+1:n_z);
        dxsum = dxsum - sum(dxkill,2);
        dx2sum = dx2sum - sum(dxkill.^2,2);
      case 'none'
        % do nothing
      otherwise
        error(['backtrack mode ' options.backtrack ' not recognized']);
    end
  end
  if options.plot_criterion
    loglog([sq_summed_disp(:) summed_sq_disp(:)])
    if (n_z < N)
      % Place a visual marker to show criterion was satisfied
      line(n_z,[sq_summed_disp(n_z) summed_sq_disp(n_z)],'LineStyle','none','Marker','x')
    end      
    yl = get(gca,'YLim');
    firstnz = find(summed_sq_disp > 0,1,'first');
    set(gca,'TickDir','out','YLim',[summed_sq_disp(firstnz) yl(2)],'XLim',[0.9 1.1*n_z]);
  end
  if (n > 0)
    dxOut = dxsum/n;
    dx2Out = dx2sum/n;
  else
    % There were no points in the neighborhod
    dxOut = zeros(d,1);
    dx2Out = zeros(d,1);
  end
end
