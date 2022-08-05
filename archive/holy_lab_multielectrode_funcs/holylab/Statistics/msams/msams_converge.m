function [y,map,n,msams_info] = msams_converge(x,y,options)
% MSAMS_CONVERGE: take mean shift steps until convergence
% This algorithm uses Minimally-Significant Adaptive Mean Shift (MSAMS),
% optionally with a variable metric, and by default leapfrog convergence
% (although this can be disabled).
%
% Syntax:
%   yf = msams_converge(x,y)
% Perform MSAMS with fixed metric.  Here, x is a d-by-N matrix of data
% points, and y is a d-by-n_landmarks matrix of landmarks.  Each landmark
% will be moved by MSAMS until it stops moving.  Note this syntax does
% not use leapfrog convergence, and thus is slow!
%
%   [yf,map,n] = msams_converge(x,y)
% This performs MSAMS with a fixed metric, but permits leapfrog
% convergence, which is much faster.  map is an index vector assigning
% the ith landmark to the landmark map(i).  n is a vector indicating, for
% each landmark, the number of data points within its local neighborhood.
%
%   [...] = msams_converge(x,y,options)
% This allows you to control MSAMS through various fields of the
% structure "options":
%   show_progress (default false): if true, prints out the number of
%     landmarks that have been processed.
%   leapfrog (default true): expect leapfrog convergence. You have to
%     request the "map" and "n" outputs for this preference to be
%     respected.
%   convergence_thresh (default eps): the tolerance for detecting cycling.
%   variable_metric (default false): if true, uses a variable metric in
%     moving landmarks (see VMAMS_STEP1). If you get warnings suggesting
%     that you use VMAMS, but don't want to, you can suppress these
%     warnings by setting options.variable_metric = false on input.
%   use_experience (default true): when true, the initial guess for the
%     scale is based on the scales found for previous landmarks (relevant
%     only when variable_metric = true).
%   debug (default false): provide graphical output showing the progress.
%     Works only in 2d.
%
% Be aware that options also gets passed to MSAMS_STEP1 or VMAMS_STEP1,
% so all of those fields apply as well.
%
%   [yf,map,n,msams_info] = msams_converge(...)
% also returns various diagnostic information.
%
% See also: MSAMS, MSAMS_STEP1, VMAMS_STEP1.
  
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
  % Test for trivial input
  if (N == 0)
    error('No data points supplied');
  elseif (N == 1)
    y = repmat(x,[1 nlandmarks]);
    map = ones(1,nlandmarks);
    n = map;
    % msams_info?
    return
  end
  
  use_weighted_metric = false;
  if isfield(options,'metric')
    use_weighted_metric = true;
    options.variable_metric = false;
  end
  if ~isfield(options,'variable_metric')
    options.variable_metric = false;
    % Check to see whether we should suggest to the user that
    % 'variable_metric' should be turned on. Note this gets skipped if
    % the user explicitly sets options.variable_metric.
    coord_var = var(x,0,2);
    min_coord_var = min(coord_var(coord_var > 0));
    max_coord_var = max(coord_var);
    if (10*min_coord_var < max_coord_var)
      warning('MSAMS:useVariableMetric',...
	      ['Large width differences among coordinates detected. Consider scaling the\n' ...
        'coordinates and/or setting variable_metric to true. To disable this\n' ...
        'warning, set variable_metric to false.'])
    end
  end
  if ~isfield(options,'leapfrog')
    options.leapfrog = (nargout > 1);
  end
  if ~isfield(options,'use_experience')
    options.use_experience = true;
  end
  if ~isfield(options,'convergence_thresh')
    % threshold for determining cycling (eps = very strict) or convergence
    if options.leapfrog
      options.convergence_thresh = eps;
    else
      options.convergence_thresh = 1e-8;
    end
  end
  if ~isfield(options,'max_iter_msams')
    options.max_iter_msams = 200;
  end
  if ~isfield(options,'show_progress')
    % If true, prints out # of moving landmarks on each iteration
    options.show_progress = false;
  end
  if ~isfield(options,'save_trajectories')
    options.save_trajectories = false;
  end
  if ~isfield(options,'save_nhistory')
    options.save_nhistory = false;
  end
  if (nargout > 3)
    msams_info = struct;
  end
% $$$   if (~isa(x,'double') || ~isa(y,'double'))
% $$$     if (isfield(options,'use_mex') && options.use_mex)
% $$$       error('MSAMS:notDouble','To use the MEX file, x and y must be of class double');
% $$$     elseif ~isfield(options,'use_mex')
% $$$       warning('MSAMS:notDouble','Not using MEX file because at least one of x or y is not of class double');
% $$$     end
% $$$     options.use_mex = false;
% $$$   end
% $$$   if ~isfield(options,'use_mex')
% $$$     options.use_mex = (exist('adaptive_meanshift_step','file') == 3);
% $$$     if ~options.use_mex
% $$$       warning('MSAMS:noMEXfile','MEX file "adaptive_meanshift_step" is not available. You might consider compiling it.')
% $$$     end
% $$$   end
% $$$   if ~isfield(options,'factor')
% $$$     % factor = alpha, the number of times by which the bias has to exceed
% $$$     % the standard error before significance is determined
% $$$     options.factor = 3;
% $$$   end
% $$$   if ~isfield(options,'backtrack')
% $$$     options.backtrack = true;
% $$$   end
  if isfield(options,'debug')
    debug = options.debug;
  else
    debug = false;
  end
  if (debug)
    figure
  end
  
  
  
  %
  % Move points under MSAMS (Minimally-Significant Adaptive Mean Shift)
  %
  y0 = y;
  map = 1:nlandmarks;
  did_not_converge = 0;
  if options.show_progress
    last_time = cputime;
  end
  for landmarkIndex = 1:nlandmarks
    if options.show_progress
      if (cputime < last_time || last_time+5 < cputime)
        last_time = cputime;
        fprintf('%d..',landmarkIndex);
      end
    end
    ycur = y0(:,landmarkIndex);
    if options.variable_metric
      if options.use_experience && isfield(msams_info,'scale')
        % Use previous iterations to get a guess for the scale.
        % This isn't necessary, but it's helpful.
        past_scales = msams_info.scale;
        past_scales_norm = sum(past_scales.^2,1);
        scale = median(past_scales ./ repmat(past_scales_norm,[d 1]),2);
      else
        scale = ones(d,1);
      end
    end
    ytraj = zeros(d,0);  % for testing cycling & reporting the history
    ntraj = [];          % for assigning correct priority in case of cycling
    scaletraj = zeros(d,0);  % for holding history of scaling in var. metric
    isdone = false;
    ismapped = false;
    isstopped = false;
    iter = 0;
    if debug
      plot(x(1,:),x(2,:),'.'); hold on; plot(ycur(1),ycur(2),'rx'); hold off; title(['Starting ' num2str(landmarkIndex)])
      pause
    end
    while ~isdone
      % Move the landmark
      ytraj(:,end+1) = ycur;
      if options.variable_metric
        % Variable metric case
        options.scale = scale;
        [ycur,neighborIndex,ncur,scale] = vmams_step1(x,ycur,options);
        % Calculate the analog of the msd output of msams_step1
        msd = scale;
        msd(msd == 0) = inf;  % to protect against divide by zero
        msd = 1./msd.^2;
        scaletraj(:,end+1) = scale;
      else
        % Fixed metric case
        if use_weighted_metric
          [sd,weights] = options.metric(ycur,x,options.metric_options);
          options.metric_weights = weights;
        end
        [ycur,neighborIndex,ncur,msd] = msams_step1(x,ycur,options);
      end
      if debug
        plot(x(1,:),x(2,:),'.'); hold on; plot(x(1,neighborIndex(1:ncur)),x(2,neighborIndex(1:ncur)),'g.'); plot(ycur(1),ycur(2),'rx'); hold off; title(num2str(landmarkIndex))
        pause
      end
      ntraj(end+1) = ncur;
      % Find the landmark starting location closest to ycur, to see
      % whether we can assign this landmark to another one
      if (options.leapfrog)
        if options.variable_metric
          dy = (repmat(ycur,1,nlandmarks)-y0).*repmat(scale,1,nlandmarks);
          %dy2 = dy.^2;
          % In the case of variable metric, we need to make sure that the
          % other landmark is within range to influence the current
          % landmark
          dx = (ycur - x(:,neighborIndex(ncur))).*scale; % most distant influence
          %dy2(dy2 > 1) = inf;
          %dist = sum(dy2,1);
          dist = sum(dy.^2,1);
          dist(dist > sum(dx.^2)) = inf;
          [md,mdIndex] = min(dist);
        else
          [md,mdIndex] = mindist(ycur,y0);
        end
      end
      % Can we assign this landmark to another one?
      if (options.leapfrog && mdIndex ~= landmarkIndex && ~isinf(md))
        if debug
          disp([ycur, y0(:,[landmarkIndex mdIndex])])
        end
        isdone = true;  % the landmark maps to a different landmark
        ismapped = true;
        map(landmarkIndex) = mdIndex;
      else
        % Can't quit yet. But maybe we've entered a cycle, in which case
        % we should quit anyway
        dy = repmat(ycur,1,size(ytraj,2)) - ytraj;
        issame = all(dy.^2 <= options.convergence_thresh*repmat(msd,[1 size(dy,2)]),1);
        if any(issame,2)
          isdone = true;
          if debug
            fprintf('Cycling...\n');
          end
          isstopped = issame(end);  % when a cycle means convergence
        end
      end
      iter = iter+1;
      isdone = isdone | iter >= options.max_iter_msams;
    end
    n(landmarkIndex) = ntraj(end);
    if (~ismapped && ~isstopped)
      % When we terminated the cycle, we might have done so at a point
      % that was not the "best" we've seen in terms of large n. Restore
      % that setting, so we stay at the best peak.
      ntrajtmp = ntraj;
      ntrajtmp(1:find(issame,1,'first')-1) = 0;  % consider history only
      % since cycling started
      [max_ntraj,max_ntraj_index] = max(ntrajtmp);
      ycur = ytraj(:,max_ntraj_index);
      n(landmarkIndex) = max_ntraj;
      if options.variable_metric
        scale = scaletraj(:,max_ntraj_index);
      end
    end
    % Do recordkeeping
    y(:,landmarkIndex) = ycur;
    if options.save_trajectories
      msams_info.ytraj{landmarkIndex} = ytraj;
    end
    if options.save_nhistory
      msams_info.nhistory{landmarkIndex} = ntraj;
    end
    if options.variable_metric
      msams_info.scale(:,landmarkIndex) = scale;
    end
    if (iter >= options.max_iter_msams)
      did_not_converge = did_not_converge+1;
    end
  end
  if options.show_progress
    fprintf('\n');
  end
  if did_not_converge
    warning('MSAMS:noConvergence','%d landmarks did not converge', ...
      did_not_converge);
  end