function [T2,T2ls] = cn_mcT2linesearch(d,nlist,params,nsim)
% cn_mcT2linesearch: Monte Carlo simulations of T^2 during linesearches
%
% Syntax:
%   [T2,T2ls] = cn_mcT2linesearch(d,nlist,params)
%   [T2,T2ls] = cn_mcT2linesearch(d,nlist,params,nsim)
% where
%   d is the dimensionality
%   nlist is a vector containing the different neighborhood sizes you want
%     to analyze, e.g., nlist = unique(round(logspace(0,3,31)));
%   params has fields "covarianceModel" and, if nsim is not supplied,
%     "pvalue". You may also supply "sim_nfactor" (with a default value of
%     2*d) to control the # of data points in the fake data sets relative to
%     max(nlist).
%   nsim controls the number of simulations run. If this is not supplied
%     explicitly, nsim = ceil(1/params.pvalue^2);
%
% The analysis time and storage requirements are O(length(nlist)*nsim), so
% be cautious about supply very low pvalues.

% Copyright 2011 by Timothy E. Holy

  %% Parse inputs
  if isscalar(nlist) || ~issorted(nlist)
    error('Must provide the list of n-values for which you want to evaluate the threshold, in sorted order');
  end
  n_keep = length(nlist);
  params = default(params,'sim_nfactor',2*d);
  if (nargin < 4)
    nsim = ceil(1/params.pvalue^2);
  end
  
  %% Prepare output storage
  T2 = inf(nsim,n_keep);
  T2ls = inf(nsim,n_keep);
  
  %% Run simulations
  for i = 1:nsim
    % Report on progress
    if (mod(i-1,round(nsim/100))==0)
      fprintf('%d%%...',round(100*(i-1)/nsim));
    end
    % Generate data points
    x = randn(d,round(max(nlist)*params.sim_nfactor));
    % Put them in order of increasing distance from zero
    [sx2,sortOrder] = sort(sum(x.^2,1));
    % Analyze neighborhoods of all chosen sizes
    for j = 1:n_keep
      thisn = nlist(j);
      if thisn == 1 || (strcmp(params.covarianceModel,'full') && thisn <= d)
        continue;
      end
      thisx = x(:,sortOrder(1:thisn));
      % Compute T^2
      [T2(i,j),mu] = mcT2ls_T2(thisx,params.covarianceModel);
      % Pick a projection direction to mimic line search
      munorm = sqrt(sum(mu.^2));
      if (munorm ~= 0)
        muhat = mu/munorm;  % direction is along the mean
      else
        % pick first coordinate as the projection direction
        muhat = zeros(d,1);
        muhat(1) = 1;
      end
      % Represent each point by components parallel and perpendicular to
      % the line; in particular, we only need the perpendicular distance,
      % not all of the perpendicular coordiantes.
      % (This can save time if d is large.)
      proj = muhat'*x;        % the projection in the chosen direction
      proj = proj(sortOrder); % put them in order of increasing distance from zero
      perp = sqrt(sx2 - proj.^2);  % the perpendicular distance
      leftflag = proj < 0;  % points in direction of -alpha
      leftflag(1:thisn) = false;  % points not in the nbrhood
      rightflag = proj > 0; % points in direction of +alpha
      rightflag(1:thisn) = false;
      yleft = [proj(leftflag); perp(leftflag)];
      yright = [proj(rightflag); perp(rightflag)];
      ycenter = [proj(1:thisn); perp(1:thisn)];
      % Find the positions along the line, in each direction, for which the
      % distance to the farthest "inside" point is equal to the distance to
      % the closest "outside" point. This defines the span of movements
      % along the line for which the inside points are the closest "thisn"
      % points.
      alphaleft = mcT2ls_alpha(ycenter,yleft,-1);
      alpharight = mcT2ls_alpha(ycenter,yright,1);
      % Compute T2 for each of these points and pick the largest
      T2left = mcT2ls_T2(bsxfun(@minus,thisx,alphaleft*muhat),params.covarianceModel);
      T2right = mcT2ls_T2(bsxfun(@minus,thisx,alpharight*muhat),params.covarianceModel);
      T2ls(i,j) = max(T2left,T2right);
    end
  end
  fprintf('done.\n');
end

function [T2,mu] = mcT2ls_T2(x,covarianceModel)
  [d,n] = size(x);
  mu = mean(x,2);
  dx = bsxfun(@minus,x,mu);
  switch covarianceModel
    case 'isotropic'
      C = sum(dx(:).^2)/(n*d-1);
      T2 = n*sum(mu.^2)/C;
    case 'diagonal'
      C = sum(dx.^2,2)/(n-1);
      T2 = n*sum(mu.^2./C);
    case 'full'
      if (n > d)
        C = (dx*dx')/(n-1);
        T2 = n*mu'*(C\mu);
      else
        T2 = inf;
      end
  end
end

function alpha = mcT2ls_alpha(ycenter,youter,s)
%   figure; plot(ycenter(1,:),ycenter(2,:),'b.'); hold on; plot(youter(1,:),youter(2,:),'r.'); plot(0,0,'kx'); hold off; axis equal;
  if (s < 0)
    [~,nbrs] = max(ycenter(1,:)); % rightmost of the inside points
  else
    [~,nbrs] = min(ycenter(1,:)); % leftmost of the inside points
  end
  tmp = find(s*youter(1,:) > s*ycenter(1,nbrs),1,'first');  % closest of the outer points that is in the correct direction
  if ~isempty(tmp)
    nbrs(2) = tmp;
  else
    alpha = s*inf;
    return
  end
  oldnbrs = [];
  while ~isequal(nbrs,oldnbrs)
    oldnbrs = nbrs;
    y1 = ycenter(:,nbrs(1));
    y2 = youter(:,nbrs(2));
    alpha = (sum(y2.^2) - sum(y1.^2))/(2*(y2(1)-y1(1)));
    ymid = [alpha; 0]; % candidate midpoint
    dy = bsxfun(@minus,ycenter,ymid);
    [~,nbrs] = max(sum(dy.^2)); % most distant of the inside points
    dy = bsxfun(@minus,youter,ymid);
    [~,nbrs(2)] = min(sum(dy.^2)); % closest of the outside points
  end
end
