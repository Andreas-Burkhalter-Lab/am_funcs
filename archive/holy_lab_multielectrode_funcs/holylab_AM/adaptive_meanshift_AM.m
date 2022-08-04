function [clust,clustinfo] = adaptive_meanshift_AM(x,y,options,display_mean_n)
% ADAPTIVE_MEANSHIFT: perform locally-adaptive meanshift (LAMS)
%%% AM edited 8/10/15 
%%% -added option to not display 'Mean n:'
% Syntax:
%   clust = adaptive_meanshift(x,y)
%   [clust,clustinfo] = adaptive_meanshift(x,y,options)
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   y is a d-by-nlandmarks matrix of landmark points in d-dimensional
%     space;
% and
%   clust is a 1-by-nlandmarks vector, giving the cluster number
%     assigned to each landmark.
%   clustinfo is an output structure with the following fields:
%     yf is a d-by-nlandmarks matrix, with each column containing the
%       final resting point of each landmark (either the point of
%       convergence, or the point at which leapfrog took charge).
%     map is a vector, map(i) is the index of the landmark that the ith
%       landmark "flows" to when flowing uphill by adaptive_meanshift steps;
%     progress is a vector containing the number of unmapped points for each
%       iteration of the adaptive_meanshift algorithm.
%     n_history, y_traj.
%
% When you are using leapfrog convergence, a facsimile of what would have
% been the final resting position of each landmark can be obtained from
%    cfinal = clustinfo.yf(:,clustinfo.map);
%
% See also: MEANSHIFT, RMEANS, ADAPTIVE_MEANSHIFT1.
  
% Old output:
%   comp is a cell array of length nclusters, each entry containing the
%     indices of the landmarks belonging to a particular cluster;

% Copyright 2006 by Timothy E. Holy
  
  % Check the inputs and supply default options
  [d,N] = size(x);
  [dy,nlandmarks] = size(y);
  if (d ~= dy)
    error(['The dimensionality of the data and the landmarks are not ' ...
           'consistent']);
  end
  if (nargin < 3)
    options = struct;
  end
  if (~isa(x,'double') || ~isa(y,'double'))
    if (isfield(options,'use_mex') && options.use_mex)
      error('MSAMS:notDouble','To use the MEX file, x and y must be of class double');
    elseif ~isfield(options,'use_mex')
      warning('MSAMS:notDouble','Not using MEX file because at least one of x or y is not of class double');
    end
    options.use_mex = false;
  end
  if ~isfield(options,'protect_n')
    options.protect_n = 0;
  end
  if ~isfield(options,'backtrack')
    options.backtrack = true;
  end
  if (options.protect_n || ~options.backtrack)
    options.use_mex = false;
  end
  if ~isfield(options,'use_mex')
    options.use_mex = (exist('adaptive_meanshift_step','file') == 3);
    if ~options.use_mex
      warning('MSAMS:noMEXfile','MEX file "adaptive_meanshift_step" is not available. You might consider compiling it.')
    end
  end
  if ~isfield(options,'use_leapfrog')
    options.use_leapfrog = true;
  end
  if ~isfield(options,'save_trajectories')
    options.save_trajectories = false;
  end
  if ~isfield(options,'factor')
    % factor = alpha, the number of times by which the bias has to exceed
    % the standard error before significance is determined
    options.factor = 3;
  end
  if ~isfield(options,'minN')
    options.minN = 10;
  end
  if ~isfield(options,'thresh')
    % threshold for determining cycling (eps = very strict) or convergence
    if options.use_leapfrog
      options.thresh = eps;
    else
      options.thresh = 1e-4;
    end
  end
  if ~isfield(options,'maxiter')
    options.maxiter = 200;
  end
  if ~isfield(options,'showprogress')
    % If true, prints out # of moving landmarks on each iteration
    options.showprogress = false;
  end
  if ~isfield(options,'ploteach')
    % Allow graphical display on each step (useful only in 2d)
    options.ploteach = false;
  end
  if ~isfield(options,'plotproject')
    % For plotting, set up a default set of projections (which simply
    % chooses the first two dimensions)
    dlim = min(2,d);   % Plotting in 1 or 2 dimensions
    options.plotproject = zeros(dlim,d);
    options.plotproject(1:dlim,1:dlim) = eye(dlim);
  end
  if ~isfield(options,'plotpause')
    % If plotting, one can pause to examine the results after each step
    options.plotpause = false;
  end
  
  
  %
  % Move points under MSAMS (Minimally-Significant Adaptive Mean Shift)
  %
  is_not_mapped = true(1,nlandmarks);
  map = 1:nlandmarks;
  is_cycling = false(1,nlandmarks);
  is_moving = is_not_mapped & ~is_cycling;
  is_movingOld = is_moving;
  R2 = zeros(1,nlandmarks);
  n = nan(1,nlandmarks);
  y0 = y;    % save the originals, for the purpose of mapping
  yOld = zeros([d,0,0],class(y));   % Convergence testing: check for cycles
  if options.save_trajectories
    y_traj = y;
  end
  iter = 0;
  progress = 0;  % progress will be sum(is_not_mapped)
  is_done = false;
  n_history = zeros([0 nlandmarks]);
  % Associate each landmark with the nearest one to the point where it
  % gets mapped
  size(x)
  while ~is_done
    if options.use_mex
      [y(:,is_moving),R2(is_moving),n(is_moving)] = ...
        adaptive_meanshift_step(x,y(:,is_moving),options.factor,options.minN);
    else
      for i = find(is_moving)
        options.nOld = n(i);
        [y(:,i),R2(i),indexContrib] = adaptive_meanshift1(x,y(:,i),options);
        n(i) = length(indexContrib);
      end
    end
    
    %% AM 8/10/15 added conditional
    if exist('display_mean_n','var') && display_mean_n
        fprintf('Mean n: %d\n',round(mean(n(is_moving))));
    end
  %%
        if options.save_trajectories
      y_traj(:,:,end+1) = y;
    end
    n_history(end+1,:) = n;
    if options.use_leapfrog
      % For each moving landmark i, assign it to another landmark j if j's
      % starting position is closer to i's current position than i's own
      % starting point
      % Note: if we're not using leapfrog, then all points remain unmapped
      [dist,map(is_moving)] = mindist(y(:,is_moving),y0);
      is_not_mappedOld = is_not_mapped;
      is_not_mapped = (map == 1:length(map)); % points to self -> not mapped
      is_moving = is_not_mapped & ~is_cycling;
    end
    is_cyclingOld = is_cycling;
    if isequal(is_moving,is_movingOld)
      % The identity of moving points did not change from the last
      % iteration. Test to see if remaining points have entered a cycle.
      y_moving = y(:,is_moving);
      if (~isempty(y_moving) && size(yOld,2) == sum(is_moving))
        % Check for equality with a previous location
        ydiff = repmat(y_moving,[1 1 size(yOld,3)]) - yOld;
        sd = sum(ydiff.^2,1);
        sd_min = min(sd,[],3);
        sdcompare = options.thresh*R2(is_moving);
        % Here is the test to see if a point is (essentially) back to one
        % of the locations on its history
        is_cycling(is_moving) = (sd_min <= sdcompare);
        yOld(:,:,end+1) = y(:,is_moving); % append to the history
      else
        % This case should only occur once, when yOld goes from being empty
        % to being initialized.
        yOld = y(:,is_moving);
      end
    elseif ~isempty(yOld)
      % The number of moving points changed. Discard any history on the
      % points that stopped moving
      yOld = yOld(:,is_moving(is_movingOld),:);
    end
    % Update the is_moving criterion (again, if using leapfrog)
    is_movingOld = is_not_mapped & ~is_cyclingOld;
    is_moving = is_not_mapped & ~is_cycling;
    iter = iter+1;
    progress(iter) = sum(is_moving);
%     if (iter > 1 && progress(iter-1) == 1)
%       line([yOld(1,:,end), y(1,is_not_mapped)],...
%         [yOld(2,:,end),y(2,is_not_mapped)]);
%       th = linspace(0,2*pi,100);
%       r = sqrt(R2(is_not_mapped));
%       line(y(1,is_not_mapped) + r*cos(th),y(2,is_not_mapped) + r*sin(th),'Color','r')
%       pause(0.1)
%     end
    if options.showprogress
      fprintf('  %d\n',progress(iter));
    end
    if options.ploteach
      meanshiftplot(x,y,sqrt(R2),options);
    end
    %is_done = ~((~isequal(is_not_mapped,is_not_mappedOld) || ...
    %  ~all(is_cycling(is_not_mapped))) && sum(is_not_mapped) > 0 && iter < options.maxiter);
    is_done = sum(is_moving) == 0 || iter >= options.maxiter;
  end
  if (iter == options.maxiter && ~all(~is_not_mapped | is_cycling))
    warning('MSAMS:earlyExit','Iterations terminated before all points were mapped or determined to be cycling; more clusters may have been found than are truly present');
  end
  if (nargout > 1)
    clustinfo.yf = y;
    clustinfo.progress = progress;
    clustinfo.n_history = n_history;
    if options.save_trajectories
      clustinfo.y_traj = y_traj;
    end
  end
  
  %
  % Begin leapfrog convergence
  %
  
  %E = [(1:length(map))' map'];  % the directed graph
  %compE = grComp(E);  % the component of each vertex
  %[comp2,ncomp2] = agglabel(compE)
  
  % OK, now every landmark "maps" to some landmark, either itself (usually
  % rare) or to another landmark.
  % We can accelerate convergence by iterating this map, so that points
  % leapfrog uphill.
  map0 = map;
  umap = unique(map);
  umapOld = [];
  while ~isequal(umapOld,umap)
    umapOld = umap;
    map = map(map);
    umap = unique(map);
  end
  % The leapfrogging process can cut off some points: consider a map of
  % only three points, where the mapping is
  %    2 3 2
  % i.e., 1->2, 2->3, 3->2. Note there is a cycle. If we iterate this map
  % we get
  %    3 2 3
  % which is stable under further iterations with itself. However, it looks
  % like there are two "clusters," [1 3] and [2]. But in fact they were all
  % connected originally. Consolidate these splits due to cycles by
  % checking the original map again, and moving each point to our best
  % guess for the most uphill point of two possible mappings. (We could
  % just choose, say, the lower index of two possible mappings, but in fact
  % the map is interesting in its own right, and we might as well try to
  % preserve as much of its character as possible.) 
  % This usually converges very quickly, although there may not be
  % guarantees that it will do so.
  is_changed = true;
  iter = 0;
  %mapOld = [];
  %while ~isequal(map,mapOld)
  while (any(is_changed) && iter < options.maxiter)
    %mapOld = map;
    %map = min([map; map(map0)]);
    map_compare = [map; map(map0)];
    is_changed = (map_compare(1,:) ~= map_compare(2,:));
    map_compare = map_compare(:,is_changed);
    n1 = n(map_compare(1,:));  % estimate max density from landmark that ...
    n2 = n(map_compare(2,:));  % ... averages over the most data points
    [tmp,index_max] = max([n1; n2]);
    is_tie = (n1 == n2);
    if any(is_tie)
      % Ties can result it an infinte loop. Break them by choosing the one
      % with the lower index
      [tmp,index_max(is_tie)] = min(map_compare(:,is_tie));
    end
    map_uphill_index = sub2ind([2 length(n1)],index_max,1:length(n1));
    map(is_changed) = map_compare(map_uphill_index);
    iter = iter+1;
  end
  if (iter >= options.maxiter)
    error('Map repair failed to be fast')
  end
  
  % Now each cluster is defined by landmarks that map to themselves; these
  % should correspond in some sense to the peaks (or at least flattest
  % regions) of density. We need to check these to make sure they're not
  % just some "local cycle," by testing to see whether a plurality of the
  % points within the local neighborhood are a member of this cluster.
  maps_to_self = find(map == 1:length(map));
  [mind,x_to_y_map] = mindist(x,y); % assignment of points to landmarks
  sd = sqrdist(x,y(:,maps_to_self));
  n_reduced = 0;
  for clusterIndex = 1:length(maps_to_self)
    clusterlabel = maps_to_self(clusterIndex);
    nbrFlag = sd(:,clusterIndex) <= R2(clusterlabel);
    nbrmap = map(x_to_y_map(nbrFlag)); % the clusterlabel for each data nbr
    [unbrmap,tmp,relabelled_nbrmap] = unique(nbrmap);
    [tmp,multiplicity] = agglabel(relabelled_nbrmap);
    [tmp,max_multiplicity_index] = max(multiplicity);
    if (unbrmap(max_multiplicity_index) ~= clusterlabel)
      % There was another cluster that had more neighbors. Change the
      % membership
      map(clusterlabel) = unbrmap(max_multiplicity_index);
      n_reduced = n_reduced+1;
    end
  end
  %n_reduced
  % Iterate to get final assignments
  mapOld = [];
  while ~isequal(map,mapOld)
    mapOld = map;
    map = map(map);
  end
  % OK, re-label the clusters 1,2,3,...
  [umap,tmp,clust] = unique(map);
  %[comp,ncomp] = agglabel(clustnum);
  %isequal(sort(ncomp),sort(ncomp2))
%   co = get(gca,'ColorOrder');
%   cla
%   for i = 1:length(clabel);
%     line(x(1,clabel{i}),x(2,clabel{i}),'Color',co(i,:),'LineStyle','none','Marker','x');
%   end
  if (nargout > 1)
    clustinfo.map = map;
  end
  return
  

  
function meanshiftplot(x,y,R,options)
  if (isscalar(R) && size(y,2) > 1)
    R = repmat(R,1,size(y,2));
  end
  xp = options.plotproject*x;
  yp = options.plotproject*y;
  dim = size(xp,1);  % # of plot dimensions
  odim = size(x,1);  % # of original dimensions
  if (dim > 1)
    plot(xp(1,:),xp(2,:),'.')
    hold on
    angle = linspace(0,2*pi,200);
    for i = 1:size(yp,2)
      xc = yp(1,i) + R(i)*cos(angle)*sqrt(dim/odim);
      yc = yp(2,i) + R(i)*sin(angle)*sqrt(dim/odim);
      plot(xc,yc,'r')
    end
    hold off
    axis equal
  else
    hist(xp,ceil(sqrt(size(x,2))));
    hold on
    scatter(yp,zeros(1,size(y,2)),'ro');
    hold off
  end
  if isfield(options,'plottitle')
    title(options.plottitle);
  end
  if options.plotpause
    pause
  else
    drawnow
  end
  
  
  % Junk
        %
      % Why we test sum(is_not_mapped) < nlandmarks: sometimes, especially if
      % landmarks come from kmeans, it takes a cycle or two to get the
      % points moving. Don't waste memory & time in that case.
