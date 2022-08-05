function [clust,clustinfo] = msams(x,varargin)
% MSAMS: clustering by Minimally-Significant Adaptive Mean Shift
% 
% MSAMS is a variant of a mean shift clustering algorithm that adaptively
% chooses the size of the region to average over. It also employs a number
% of techniques to increase its efficiency.
%
% Syntax:
%
%   clust = msams(x)
%
% Here x is a d-by-N matrix of data points, where each column of x is one
% point in d dimensions (N points total). On output, clust will be a
% 1-by-N vector giving the cluster # assigned to each point.  If x has many
% points, you will want to read below about performance.
%
%
%   clust = msams(x,y)
%
% An alternative syntax in which a set of "probe points" y flow on a fixed
% background of "data points" x. This is particularly useful in
% bootstrap/cross-validation type studies (see below). If you do not supply
% y, by default each data point is used as a probe point.
%
%
%   clust = msams(...,lminfo)
%
% Use landmarks to help speed up the neighbor-finding operations. lminfo is
% a structure of the type returned by CHOOSE_LANDMARKS. As a rough
% guideline, you might want to choose the # of landmarks as sqrt(N), which
% for most probe points should convert the neighbor-finding algorithm to
% something like O(sqrt(N)) rather than O(N).
% By default, an lminfo is calculated for you if you do not supply it.
% (If desired, this can be turned off with the option
%     'compute_lminfo_if_missing' = false,
% but other than for debugging purposes there are presumably few reasons to
% turn it off.) Consequently, this syntax is useful mainly if you have
% already gone to the trouble to define landmarks and don't want that work
% re-done for you.
%
%
%   [clust,clustinfo] = msams(...)
%
% This provides extra diagnostic information in the structure clustinfo.
% The precise form of this structure depends on the version of the
% algorithm being run (MEX or non-MEX). The most detailed description can
% be found in the help for MSAMS_CONVERGE_MEX.
% Of particular note are the fields named "map" or "mapN" (N an integer),
% which specify the index of the data point (in x) closest to the final
% flow position of each probe point. "map" is the "final" version, and
% "mapN" represents this map at various stages of analysis. See the code
% for complete details.
%
%
%   clust = msams(...,options)
%
% This syntax allows you to control the algorithm with some options, of
% which the most important are:
%   flow_only_landmarks (performance-related, default false): Only relevant
%     when y is not supplied explicitly. If false, each data point is used
%     as a probe point. If true, only the landmark points are used as probe
%     points, and the clustering assignment for the rest of the points is
%     copied from the associated landmark point. Setting to true permits a
%     substantial speedup, at the possible risk of lower accuracy. Note
%     that when using this option, for convenience the clusters are
%     reported for the data points rather than the landmarks-as-y.
%   leapfrog (performance-related, default true): Allows early termination
%     of flow, which can greatly speed the algorithm. Applicable only if
%     either all data points are being used as probe points, or if all
%     landmarks are being used as probe points.
%     leapfrog is used to calculate a more specific setting, 'terminate_mode'
%     ('x' = terminate on changed nearest data point, 'l' = terminate on
%     changed nearest landmark, 'c' = terminate only when cycling). If you
%     prefer, you can specify this directly, in which case it takes
%     precedence.
%   reflow (default false): If flow terminates early (see leapfrog), there
%     are circumstances in which the early termination prevents accurate
%     determination of the true "final" destination of each probe point. If
%     reflow==true, each "self-mapping" probe point is re-subjected to
%     msams, this time without early termination. Any clusters whose probe
%     points move to the same place are merged.
%     Setting reflow=true tends to result in fewer clusters. The choice to
%     leave the default as "false" is largely to preserve backward
%     compatibility; from a conceptual standpoint it seems that the most
%     self-consistent choice is reflow=true.
%
% Some examples:
% Basic clustering, focusing on issues affecting performance:
%   If the number N of data points is modest (say, 10^4 or fewer), there is
%   probably no reason to do anything other than
%      [clust,clustinfo] = msams(x);
%   (By default this will still use landmarking internally to speed up the
%   neighbor-finding computations, but you won't need to know about it.)
%
%   If x contains a lot of points and you can't afford the time it takes to
%   flow each data point, then first try this:
%      [clust,clustinfo] = msams(x,struct('flow_only_landmarks',true));
%   or, if you prefer to handle the landmarks explicitly,
%      n_landmarks = ceil(sqrt(N));  % a guideline, not a rule
%      lminfo = choose_landmarks(x,n_landmarks);
%      [clust,clustinfo] = msams(x,lminfo,struct('flow_only_landmarks',true));
%
%   The limitation of the above is that each point in x is assigned to the
%   cluster of its corresponding landmark, and it is possible that the
%   assignments could be too "coarse." To improve the "resolution," you
%   could increase the number of landmarks, and in many cases this might be
%   the easiest option. However, since neighbor-finding performance should
%   be optimal near n_landmarks = sqrt(N), in some circumstances increasing
%   the number of landmarks may not be desirable. An alternative is to
%   supply "extra" y, but to permit leapfrog convergence by including all
%   of the landmarks among the y:
%      n_landmarks = ceil(sqrt(N));
%      min_y = min(N,5000);
%      if (min_y > n_landmarks)
%        lminfo = choose_landmarks(x,n_landmarks);
%        % We'll choose some extra probe points, making sure to include the
%        % landmarks among them
%        yinfo = choose_landmarks(x,min_y,lminfo);
%        [clust,clustinfo] = msams(x,lminfo,yinfo.landmarks,struct('terminate_mode','l'));
%        % Assign the x to clusters according to their closest y
%        clust = clust(yinfo.landmarkAssignment);
%      else
%        [clust,clustinfo] = msams(x,struct('flow_only_landmarks',true));
%      end
%
% See also: METRIC_BIGGEST_DIMENSIONS, MSAMS_STEP1, MSAMS_CONVERGE_MEX, CHOOSE_LANDMARKS, MSAMS_RESAMPLE, MSAMS_CONVERGE, VMAMS_STEP1.

% One of the things you can do with these options is use a non-euclidean
% metric. Currently, what is supported are weighted metrics, as in
% METRIC_BIGGEST_DIMENSIONS.  You have to specify the two fields,
% 'metric' and 'metric_options' to use any metric.
%

% Copyright 2006-2011 by Timothy E. Holy
  

% Cases:
% 1. y = x: leapfrog on, terminate when closest x is different from
%    starting x. We can still use landmarks productively to triage
%    neighborhoods.
% 2. y = landmarks ('flow_only_landmarks'): leapfrog on, terminate when
%    closest datapoint is assigned to a different landmark, or when
%    closest landmark is different from self. (These two should be nearly
%    identical in terms of # of iterations.)
% 3. User-supplied y. Is it worth checking to see whether the y include the
%    landmarks? (Only if user supplies landmarks.) If this is found to be
%    true, this could be handled like #2. If the y do not include the
%    landmarks, we want to turn leapfrog off because we don't know where
%    the landmarks go.
% So really there are 3 cases:
%   a) If we know where the data points will flow, terminate once we move
%      closer to a different data point
%   b) If not a, then see whether we know where the landmarks will flow. If
%      so, terminate once we move closer to a different landmark
%   c) If neither, then don't terminate flow early.
% But we also need to let the user turn all this off. So reserve 'leapfrog'
% for the user's control over the process, but use it to set another
% variable that decides termination: x, landmark, or none.

  %% Parse the inputs
  N = size(x,2);
  options = struct;
  have_lminfo = false;
  user_supplied_y = false;
  
  for curarg = 1:length(varargin)
    if isnumeric(varargin{curarg})
      % user-supplied y
      y = varargin{curarg};
      user_supplied_y = true;        
    elseif isstruct(varargin{curarg})
      if isfield(varargin{curarg},'landmarkAssignment')
        % this looks like a landmark structure
        lminfo = varargin{curarg};
        have_lminfo = true;
      else
        % a structure not recognized, so assumed to be options
        options = varargin{curarg};
      end
    end
  end
  
  %% Set main options
  options = default(options,'compute_lminfo_if_missing',true,...
    'flow_only_landmarks',false,...
    'leapfrog',true,...
    'consolidate',true,...
    'reflow',false);

  %% Generate landmarks if needed
  if ~have_lminfo && options.compute_lminfo_if_missing
    % One might be tempted to see whether we can use y (or some of the y)
    % for the landmarks, but there is no guarantee that these are
    % well-distributed in space. So generate them from scratch.
    lminfo = choose_landmarks(x,ceil(sqrt(N)));
    have_lminfo = true;
  end
  
  %% Determine the appropriate termination mode
  % For experts: if user_supplied_y is true, then we have to do some
  % calculations to determine whether the data points or landmarks are
  % among them. You can short-circuit this calculation by supplying a field
  % 'terminate_mode' in options (set to 'x' if all x are in y, 'l' if all
  % landmarks are in y).
  if user_supplied_y
    if options.flow_only_landmarks
      error('Supplying y is inconsistent with flow_only_landmarks = true');
    end
    if isfield(options,'terminate_mode')
      terminate_mode = options.terminate_mode;
    else
      terminate_mode = 'c';
      if options.leapfrog
        if (size(y,2) >= N)
          % See whether the probe points include all data points
          [tf,loc] = ismember(x',y','rows');
          if all(tf)
            terminate_mode = 'x';
          end
        elseif (have_lminfo && size(y,2) >= size(lminfo.landmarks,2))
          % See whether the probe points include all landmarks
          [tf,loc] = ismember(lminfo.landmarks',y','rows');
          if all(tf)
            terminate_mode = 'l';
          end
        end
      end
    end
  else
    % User did not supply y. Determine whether we are flowing all data
    % points, or just the landmarks.
    if options.flow_only_landmarks && have_lminfo
      y = lminfo.landmarks;
      if options.leapfrog
        terminate_mode = 'l';
      end
    else
      y = x;
      if options.leapfrog
        terminate_mode = 'x';
      end
    end
    if isfield(options,'terminate_mode')
      terminate_mode = options.terminate_mode; % user override
    end
  end
  options.terminate_mode = terminate_mode;  % now terminate_mode supercedes leapfrog

  %% Set other defaults
  options = default(options,'n_threads',8);
  if ~isfield(options,'use_mex')
    options.use_mex = (exist('msams_converge_mex','file') == 3) && have_lminfo;
    if isfield(options,'metric')
      options.use_mex = false;
    end
  end
  if options.reflow && ~options.use_mex
    error('Currently, reflow is supported only if using the MEX version')
  end
  

  %% Do the MSAMS steps
  if options.use_mex
    clustinfo = msams_converge_mex(x,lminfo,y,options);
    yf = clustinfo.yf;
    map = clustinfo.closestDataIndex;
    n = clustinfo.n;
  else
    [yf,map,n,clustinfo] = msams_converge(x,y,options);
    clustinfo.yf = yf;
  end  
  clustinfo.y0 = y;
  clustinfo.map0 = map;
  clustinfo.n = n;
  % Wait to set variable_metric here so that msams_converge can handle
  % defaults correctly
  options.variable_metric = isfield(clustinfo,'scale');
    
  %% In case of early termination, aggregate points
  % Three steps:
  %   1. Leapfrog flow
  %   2. "Consolidation" by voting
  %   3. Reflow
  % The indexing gets a little tricky in cases where we are not flowing all
  % of the data points. In trying to work out how to do this, I generated a
  % list of useful indexing combinations (not certain all these are
  % right...)
  %   a) if flowing landmarks and no other points, then
  %         mapx = map(lminfo.landmarkAssignment)
  %      is a "complete" flow map among the data points, and
  %         mapl = lminfo.landmarkAssignment(map)
  %      is a "complete" flow map among the landmarks. If mapl is iterated and
  %      finalized, then
  %         mapx2x = clustinfo.map0(mapl(lminfo.landmarkAssignment))
  %      is a "final" mapping of each x. The yf of the "stable" points can
  %      be extracted as
  %         yfstable = yf(:,mapl == 1:length(mapl));
  %   b) if flowing y, then
  %      i) if y includes x, then
  %           mapy = loc(map)
  %         is a flow map among y, and
  %           mapx = map(loc)
  %         is a flow map among x. If mapx is finalized,
  %           mapy2x = mapx(clusterinfo.map0)))
  %         is a "final" mapping of each y to x. yf for the "stable" points
  %         is extracted as
  %            yf = yf(:,loc);    % just the final location for the xs
  %            yfstable = yf(:,mapx == 1:length(mapx));
  %     ii) if y includes the landmarks, then
  %           mapy = loc(lminfo.landmarkAssignment(map))
  %         is a flow map among y, and
  %           mapl = lminfo.landmarkAssignment(map(loc))
  %         is a flow map among the landmarks. If mapl is finalized,
  %           mapy2l = mapl(lminfo.landmarkAssignment(clustinfo.map0)))
  %         is the "final" mapping of y to landmarks.
  %           mapy2x = clustinfo.map0(mapy2l)
  %         is a "final" mapping of each y to x. yf for the "stable" points
  %         is extracted as
  %           yf = yf(:,loc)
  %           yfstable = yf(:,mapl == 1:length(mapl));
  %    iii) if y does not include all of the landmarks, then "cluster
  %         number" is inferred from clustinfo.map0 (note in this case
  %         terminate_mode == c)
  %     
  if (options.terminate_mode ~= 'c')
    % Do the nasty indexing gymnastics to get relationships among x, y, and
    % the "complete" map (whether that is landmarks or x)
    complete2x = clustinfo.map0;
    x2complete = 1:size(x,2);
    if (options.terminate_mode == 'l')
      % "complete map" will be among landmarks
      map = lminfo.landmarkAssignment(map);
      x2complete = lminfo.landmarkAssignment(x2complete);
    end
    y2complete = map;  % convert y indices to whatever terms the "complete map" will be in
    if user_supplied_y
      map = map(loc);  % select just among the "complete map"
      complete2x = complete2x(loc);
      n = n(loc);
      yf = yf(:,loc);
    end
    % Do leapfrog flow (first stage of finalizing the map)
    map0 = map;
    map = leapfrog(map,n);
    % Store nearest x, indexed by y, after first stage of finalizing
    clustinfo.map1 = complete2x(map(y2complete));
  
    %
    % Now each cluster is defined by points that map to themselves; these
    % should correspond in some sense to the peaks (or at least flattest
    % regions) of density. We need to check these to make sure they're not
    % just some "local cycle" by testing to see whether a plurality of the
    % points within the local neighborhood are a member of this cluster.
    % Essentially, each point identified as belonging to the local
    % neighborhood of a self-mapping landmark will "vote" for the assignment
    % of the cluster overall.
    %
    if options.consolidate
      mapOld = [];
      clustinfo.n_consolidated = 0;
      % For all that start out mapping to self, determine the list of
      % neighbors. We only need to do this once.
      maps_to_self0 = find(map == 1:length(map));
      maps_to_self = maps_to_self0;
      neighbor_list = cell(1,length(maps_to_self));
      for clusterIndex = 1:length(maps_to_self)
        tLandmark = maps_to_self(clusterIndex);
        sd = sqrdist(x,yf(:,tLandmark));
        [sd,neighbor_list_tmp] = sort(sd);
        neighbor_list{clusterIndex} = neighbor_list_tmp(1:n(tLandmark));
      end
      while ~isequal(mapOld,map)
        mapOld = map;
        n_maps_to_self = length(maps_to_self);
        votes = zeros(n_maps_to_self,n_maps_to_self);
        index0 = findainb(maps_to_self,maps_to_self0);
        for clusterIndex = 1:length(maps_to_self)
          % Recreate the local set of points, as the list "indexContrib"
          % We have to recalculate this (rather than saving it from before),
          % because cycling might mean that the last iteration was not the
          % one kept, and it's too memory intensive to store the full list on
          % each iteration.
          indexContrib = neighbor_list{index0(clusterIndex)};
          % Of the points in the neighborhood, what is their cluster
          % assignment?
          if terminate_mode == 'x'
            nbrmap = map(indexContrib);
          elseif terminate_mode == 'l'
            nbrmap = map(lminfo.landmarkAssignment(indexContrib));
          end
          % Calculate the number of points assigned to each cluster within the
          % local neighborhood---this tallies the "votes"
          [unbrmap,tmp,relabelled_nbrmap] = unique(nbrmap);
          [tmp,votes_tmp] = agglabel(relabelled_nbrmap);
          votes(clusterIndex,findainb(unbrmap,maps_to_self)) = votes_tmp;
        end
        % Re-assign each self-mapping point to the one that receives the
        % most "votes" in its local neighborhood. Note we do this only after
        % all votes have been tallied for all "candidates," so that there can't
        % be any order effects.  Not that those seem likely...
        [max_votes,maxIndex] = max(votes,[],2);
        map(maps_to_self) = maps_to_self(maxIndex);
        clustinfo.n_consolidated = clustinfo.n_consolidated + sum(maxIndex' ~= 1:n_maps_to_self);
        % Iterate map to get final assignments
        map = leapfrog(map);
        maps_to_self = find(map == 1:length(map));
      end
    end
    % Store after 2nd stage of finalizing
    clustinfo.map2 = complete2x(map(y2complete));
    
    if options.reflow
      % Restart the flow at the "peak" of each self-mapping point,
      % but without early termination
      ops = options;
      ops.terminate_mode = 'c';
      maps_to_selfOld = [];
      maps_to_self = find(map == 1:length(map));
      while ~isequal(maps_to_selfOld,maps_to_self)
        newy = yf(:,maps_to_self);
        newclustinfo = msams_converge_mex(x,lminfo,newy,ops);
        % We "label" each final flow point by the x point it is closest to.
        % This will define the unique clusters.
        stable2x = newclustinfo.closestDataIndex;
        % Now fix the "base" map with this new information
        map0(maps_to_self) = x2complete(stable2x);
        n(maps_to_self) = newclustinfo.n;
        map = leapfrog(map0,n);
        maps_to_selfOld = maps_to_self;
        maps_to_self = find(map == 1:length(map));
      end
      complete2x = map0;    
%       % Now check to see whether the data points closest to these
%       % probe-flow termini also terminate in the same place (meaning, end
%       % up being closest to themselves out of all data points)
%       [ustablex,~,uIndex] = unique(stable2x);
%       newclustinfo2 = msams_converge_mex(x,lminfo,x(:,ustablex),ops);
%       stable2xnew = newclustinfo2.closestDataIndex;
%       stable2x = stable2xnew(uIndex);
%       % Convert to complete self-mapping set
%       stable2complete = stable2x;
%       if (options.terminate_mode == 'l')
%         stable2complete = lminfo.landmarkAssignment(stable2x);
%       end
%       % Define the complete map to account for the new destinations
%       newmap = zeros(size(map));
%       newmap(maps_to_self) = stable2complete;
%       ustable2complete = unique(stable2complete);
%       newmap(ustable2complete) = ustable2complete;
%       map(ustable2complete) = ustable2complete;
%       map = newmap(map);
%       % Use the new information about truly-final flow targets to update
%       % complete2x
%       complete2x(map(maps_to_self)) = stable2x;
    end
    
    map = complete2x(map(y2complete));
  end

  clustinfo.map = map;
  % Re-label the clusters 1,2,3,...
  [umap,tmp,clust] = unique(map);
  if options.flow_only_landmarks
    clust = clust(lminfo.landmarkAssignment);
  end
