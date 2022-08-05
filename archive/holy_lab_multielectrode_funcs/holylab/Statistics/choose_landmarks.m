function lminfo = choose_landmarks(x,k,varargin)
% choose_landmarks: split a group of points up into local regions
% Syntax:
%   lminfo = choose_landmarks(x,k)
%   lminfo = choose_landmarks(x,k,lminfoOld)
%   lminfo = choose_landmarks(...,options)
% where
%   x is a d-by-N matrix of data points, each column a separate data
%     point in d dimensions;
%   k is the number of landmarks desired;
%   lminfoOld (optional) is the output of a previous call to
%     choose_landmarks; use this syntax if you want to add some new
%     landmarks (k must be larger than the # of old landmarks).
%   options (optional) is a structure with the following fields:
%     frac_random (default 0): the fraction of landmarks chosen randomly
%       from the data points. The remainder of landmarks will be chosen
%       iteratively as the point farthest from the existing landmarks.
%       A value of 0 will choose the first landmark randomly, and the rest
%       according to distance. frac_random = 0 will tend to make landmark
%       groupings of very similar physical size; frac_random = 1 will tend
%       to make landmark groupings containing relatively similar #s of
%       points.
%     metric (default @sqrdist): a handle to a function that computes the
%       distance between points.  This function must be able to handle a
%       matrix of data points for its first argument. 'sqrdist' is an
%       example of such a function.
%     seed_landmarks: if present, starts the process of landmark creation
%       using the supplied landmarks, and then creates the rest as needed.
%     center: if true, causes the landmarks to move inward from the "edge"
%       of their assigned group towards the mean (it's like doing a single
%       round of k-means).
%     progress (default true): if true, will provide an update about
%       progress every five seconds.
% and
%   lminfo is an output structure with the following fields:
%     landmarks: a d-by-k matrix, giving the positions of the landmarks;
%     landmarkAssignment: a 1-by-N vector, landmarkAssignment(i) is the
%       index of the landmark closest to the ith data point, x(:,i);
%     landmarkList: a 1-by-k cell array, each element giving the indices of
%       the points in x closest to this landmark;
%     dist_to_closest_landmark: an N-by-1 vector, containing the distance
%       to the closest landmark for each data point;
%     radius_of_neighborhood: a k-by-1 vector, giving the distance to the
%       farthest point assigned to each landmark;
%     landmark_xIndex: a 1-by-k vector, giving the indices of the points in
%       x used as landmarks (NaNs are used for seed_landmarks);
%     metric: the function handle used to compute distances.
%
% See also: REINDEX_LANDMARKS.
  
% Copyright 2006-2011 by Timothy E. Holy.
% 2011: added ability to pass in old landmark structure (TEH)

  %% Parse inputs
  options = struct;
  have_some_landmarks = false;
  for i = 1:length(varargin)
    if ~isstruct(varargin{i})
      error('Final inputs must be structures');
    end
    if isfield(varargin{i},'landmarkAssignment')
      lminfo = varargin{i};
      have_some_landmarks = true;
    else
      options = varargin{i};
    end
  end
  options = default(options,'metric',@sqrdist,...
    'frac_random',0,...
    'center',false,...
    'progress',true);
  if have_some_landmarks
    if ~isequal(options.metric,lminfo.metric)
      error('Cannot ask for new landmarks with a different metric');
    end
  end

  % See if we're using sqrdist (that will allow us to use mindist to save memory)
  funcstr = func2str(options.metric);
  if isequal(funcstr,'sqrdist')
    using_sqrdist = true;
  end
  
  [d,N] = size(x);
  
  %% Keep track of time spent
  if options.progress
    first_update = true;
    % Start the timer, but try to avoid resetting it if possible
    try
      tprev = toc;
    catch
      tic;
      tprev = 0;
    end
  end
  
  %% Handle the case when the number of landmarks is equal to data size
  if(k>=N)
    if have_some_landmarks
      if (size(lminfo.landmarks,2) == N)
        return;
      else
        newIndex = setdiff(1:N,lminfo.landmark_xIndex);
        allIndex = [lminfo.landmark_xIndex newIndex];
        lminfo.landmarks(:,end+1:N) = x(:,newIndex);
        lminfo.landmarkAssignment(allIndex) = 1:N;
        lminfo.landmarkList = num2cell(allIndex);
        lminfo.dist_to_closest_landmark = zeros(1,N);
        lminfo.radius_of_neighborhood = zeros(1,N);
        lminfo.landmark_xIndex = allIndex;
        return;
      end
    end
    lminfo.landmarks=x;
    lminfo.landmarkAssignment=1:N;
    lminfo.landmarkList=agglabel(1:N);
    lminfo.dist_to_closest_landmark = zeros(1,N);
    lminfo.radius_of_neighborhood = zeros(1,N);
    lminfo.landmark_xIndex = 1:N;
    lminfo.metric = options.metric;
    return
  end
  if (have_some_landmarks && k < size(lminfo.landmarks,2))
    error('Cannot ask for fewer landmarks than you already possess.')
  end
  
  
  %% Initialize
  if have_some_landmarks
    % The user supplied old landmarks, start with these
    n_start = size(lminfo.landmarks,2);
    y(:,1:n_start) = lminfo.landmarks;
    xIndex(1:n_start) = lminfo.landmark_xIndex;
    R2x = lminfo.dist_to_closest_landmark.^2;
    landmarkIndex = lminfo.landmarkAssignment;
    % Pre-allocate rest (for speed)
    y(:,n_start+1:k) = nan;
    xIndex(n_start+1:k) = nan;
  else
    if isfield(options,'seed_landmarks')
      % Use the seed points as the first landmarks
      n_start = min(k,size(options.seed_landmarks,2));
      y = options.seed_landmarks(:,1:n_start);
      xIndex = nan(1,n_start);
    else
      % Initialize with random point(s)
      n_start = ceil(options.frac_random * k);
      if (n_start < 1)
        n_start = 1;
      end
      if (n_start > k)
        n_start = k;
      end
      if (n_start == 1)
        xIndex = round(N*rand(1,1)+0.5);
      else
        xIndex = randperm(N);
      end
      xIndex = xIndex(1:n_start);
      y = x(:,xIndex);
    end
    % Calculate distances to "starter" landmarks
    % Do this in a way that will work even if memory is tight
    if using_sqrdist
      if (options.progress && d*N*n_start > 1e7)
        fprintf('Choosing %d landmarks: ',k);
        first_update = false;
      end
      [R2x,landmarkIndex] = mindist(x,y);
    else
      try
        cdist = options.metric(x,y);
        % For each data point, calculate the distance to the closest landmark (&
        % keep track of which one it is)
        [R2x,landmarkIndex] = min(cdist,[],2);
        landmarkIndex = landmarkIndex';
      catch
        % We go here if memory is insufficient for doing them all at once
        R2x = zeros(1,N);
        landmarkIndex = zeros(1,N);
        for i = 1:N
          cdist = options.metric(x(:,i),y);
          [R2x(i),landmarkIndex(i)] = min(cdist);
        end
      end
    end
    % In case the data set contains duplicated points, we need to make sure
    % that every one of the starter landmarks has points assigned, and toss
    % any that have no points
    [xlabel,nlabel] = agglabel(landmarkIndex);
    keepFlag = (nlabel > 0);
    xlabel = xlabel(keepFlag);
    n_start = sum(keepFlag);
    xIndex = xIndex(keepFlag);
    y = y(:,keepFlag);
    for i = 1:n_start
      landmarkIndex(xlabel{i}) = i;
    end
    % Pre-allocate remainder
    xIndex(n_start+1:k) = nan;
    y(:,n_start+1:k) = nan;
  end
  
  
  %% Choose remaining landmarks iteratively
  % Find the point that is farthest from the current landmarks
  for yIndex = n_start+1:k
    % Find the point farthest from the existing landmarks
    [mxR2,maxIndex] = max(R2x);
    if (mxR2 == 0)
      fprintf('\n');
      warning('Getting fewer landmarks, because all possibilities are already taken');
      k = yIndex-1;
      break;
    end
    xIndex(yIndex) = maxIndex;  % this is the new landmark
    % Compute the distance of each point to this new landmark
    cdist = options.metric(x(:,maxIndex),x);
    % Determine the set of points for which this landmark is the closest yet
    is_closest = (cdist < R2x);
    landmarkIndex(is_closest) = yIndex;
    R2x(is_closest) = cdist(is_closest);
    y(:,yIndex) = x(:,maxIndex);
    if options.progress
      if (toc - tprev > 5)
        tprev = toc;
        if first_update
          fprintf('Choosing %d landmarks: ',k);
          first_update = false;
        end
        fprintf('%d...',yIndex);
      end
    end
  end

  if (options.progress && ~first_update)
    fprintf('done\n');
  end

  xlabel = agglabel(landmarkIndex);
  %% Perform centering, if desired
  if (options.center && using_sqrdist)
    % Move each landmark to the center of the points that are currently
    % assigned to it
    for yIndex = 1:k
      y(:,yIndex) = mean(x(:,xlabel{yIndex}),2);
    end
    [R2x,landmarkIndex] = mindist(x,y);
    [xlabel,nlabel] = agglabel(landmarkIndex);
    % Toss any that have no points assigned
    keepflag = (nlabel > 0);
    y = y(:,keepflag);
    xlabel = xlabel(keepflag);
    k = sum(keepflag);
  end

  %% For each landmark, calculate distance to farthest assigned point
  R2y = nan(1,k);
  for yIndex = 1:k
    R2y(yIndex) = max(R2x(xlabel{yIndex}));
  end
  
  %% Put results in output form
  lminfo.landmarks = y;
  lminfo.landmarkAssignment = landmarkIndex;
  lminfo.landmarkList = xlabel;
  lminfo.dist_to_closest_landmark = sqrt(R2x);
  lminfo.radius_of_neighborhood = sqrt(R2y(:));
  lminfo.landmark_xIndex = xIndex;
  lminfo.metric = options.metric;
