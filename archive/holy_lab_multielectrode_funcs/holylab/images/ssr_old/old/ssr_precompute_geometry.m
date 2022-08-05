function features = ssr_initialize_features(input_data,previous_features, ...
					    params)
            
  % Note a general indexing convention: the order is
  %    coordinate, feature, stack#/time
  % Not all things need all these indices, of course, but they preserve
  % ordering of the relevant indices

  if ~isfield(params,'register_each_stack')
    params.register_each_stack = true;
  end
  
  filemode = iscell(input_data);
  if filemode
    % Input is a cell array of filenames
    % Parse names to get stack #s
    features.stacknum = nan(size(input_data));
    for i = 1:numel(input_data)
      [tmp,count] = sscanf(input_data{i},'stack%d.mat',1);
      if ~isscalar(tmp)
        error('Error parsing stack numbers');
      end
      features.stacknum(i) = tmp;
    end
    % Parse directory name to get pyramid structure
    s = pwd;
    ptoken = 'pyramids';
    pyramidIndex = findstr(ptoken,s);
    if ~isempty(pyramidIndex)
      if ~isscalar(pyramidIndex)
        error('Error parsing directory name for pyramid level');
      end
      s = s(pyramidIndex+length(ptoken):end);  % discard all but dimension info
      features.pyramid_dimensions = sscanf(s,'%dx')';
      if numel(features.pyramid_dimensions) ~= 3
        error('Error parsing pyramid dimensions');
      end
    end
    % Get the stack dimensions, and load that first image
    load(input_data{1});
    features.stacksize = size(im);
    features.spacing = spacing;
%     if length(features.stacksize == 2)
%       % Check to see if this is really 1d
%       features.stacksize(features.stacksize == 1) = [];
%     end
  else
    % Input is a stackmm/stack#s structure
    features.stacknum = input_data.stacknum;
    features.pyramid_dimensions = [1 1 1];
    features.stacksize = input_data.smm.size;
    h = input_data.header;
    features.spacing = [h.um_per_pixel_xy([1 1]) ...
			diff(h.piezo_start_stop)/h/frames_per_stack];
    im = input_data.smm(:,:,:,features.stacknum(1));
  end
  stack_dims = length(features.stacksize);
  n_stacks = length(features.stacknum);

  % Subdivide the space of the stack into regions; each feature will be
  % supported on one region
  features.center = params.center;
  features.margin = params.margin;
  for i = 1:stack_dims
    xmin{i} = (0:params.center(i):features.stacksize(i)-1) - params.margin(i);
    xmax{i} = xmin{i} + params.center(i) - 1 + 2*params.margin(i);
  end
  features.Xmin = grid2col(xmin)';
  features.Xmax = grid2col(xmax)';
  features.Xcenter = (features.Xmin + features.Xmax)/2;
  n_features = size(features.Xmin,2);

  % Precalculate the list of pairwise overlaps
  % This will be needed for updating T & the registration parameters
  intersect_region = cell(1,stack_dims);
  % Keep track of pairs that overlap.
  % Start with a list n_features in length, and grow in chunks as necessary
  % (more efficient than re-allocating & copying each time)
  overlappingPairsIndex = nan(2,n_features);
  overlappingPairsCoords = cell(2,n_features);
  curIndex = 1;
  for featureIndex1 = 1:n_features
    % Start featureIndex2 at featureIndex1+1 because of symmetry, and we already know that each
    % one overlaps with itself
    for featureIndex2 = featureIndex1+1:n_features
      % Check to see if support regions overlap
      overlap = true;
      dimIndex = 1;
      while (dimIndex <= stack_dims && overlap)
        intersect_region{dimIndex} = IntersectIntervals(...
          [features.Xmin(dimIndex,featureIndex1) features.Xmax(dimIndex,featureIndex1)],...
          [features.Xmin(dimIndex,featureIndex2) features.Xmax(dimIndex,featureIndex2)]);
        if isempty(intersect_region{dimIndex})
          overlap = false;
        end
        dimIndex = dimIndex+1;
      end
      if overlap
        if (curIndex > size(overlappingPairsIndex,2))
          % Grow a chunk
          overlappingPairsIndex(:,end+1:end+n_features) = nan;
          overlappingPairsCoords{2,end+n_features} = [];
        end
        overlappingPairsIndex(:,curIndex) = [featureIndex1;featureIndex2];
        for dimIndex = 1:stack_dims
          overlappingPairsCoords{1,curIndex}(dimIndex,:) = intersect_region{dimIndex}' - ...
            features.Xmin(dimIndex,featureIndex1)+1;
          overlappingPairsCoords{2,curIndex}(dimIndex,:) = intersect_region{dimIndex}' - ...
            features.Xmin(dimIndex,featureIndex2)+1;
        end
        curIndex = curIndex+1;
      end
    end
  end
  features.overlappingPairsIndex = ...
    overlappingPairsIndex(:,1:curIndex-1)';
  features.overlappingPairsCoords = ...
    overlappingPairsCoords(:,1:curIndex-1)';
  
  % Precalculate pixel blocks that have common membership in features
  % (i.e., groups of pixels that all belong to the same set of features)
  % This will be needed for updating S
  % These pixel blocks "start" at the beginning of each feature, and
  % "end" just before the beginning of the "next" feature (this works
  % cleanly in 1-d, we have to grid it in higher dimensions)
  for i = 1:stack_dims
    xstart = unique([xmin{i} xmax{i}+1]); % The last one is only an end-marker
    xblock_min{i} = xstart(1:end-1);
    xblock_max{i} = xstart(2:end)-1;
  end
  Xblock_min = grid2col(xblock_min);   % contains gridded block starts
  Xblock_max = grid2col(xblock_max);   % contains block ends
  n_blocks = size(Xblock_min,1);
  % Preallocate
  block(n_blocks).image_coords_min = [];
  block(n_blocks).image_coords_max = [];
  block(n_blocks).index = [];
  block(n_blocks).features_coords_min = [];
  block(n_blocks).features_coords_max = [];
  for i = 1:n_blocks
    % Calc the offset of the current pixel block from the beginning of
    % each feature
    min_offset = repmat(Xblock_min(i,:)',1,n_features) - features.Xmin;
    % Determine which features contain this block
    contains_block = ...
      all(min_offset >= 0,1) & ...
      all(features.Xmax>=repmat(Xblock_max(i,:)',1,n_features),1);
    featuresIndex = find(contains_block);
    % Save the index of features containing this block, as well as the
    % relative offsets of this block in each feature
    block(i).image_coords_min = Xblock_min(i,:)';
    block(i).image_coords_max = Xblock_max(i,:)';
    block(i).index = featuresIndex;
    block(i).features_coords_min = min_offset(:,featuresIndex)+1;
    block(i).features_coords_max = ...
      repmat(Xblock_max(i,:)',1,length(featuresIndex)) - ...
      features.Xmin(:,featuresIndex) + 1;
  end
  features.block = block;
  
  % Initialize registration parameters, using previous_features
  if ~isempty(previous_features)
    % Copy the old registration information
    registration = previous_features.registration;
    % Since the number of features might (is expected to) change, we have
    % to assign each old feature to a new feature.
    % Find the closest correspondence between features
    % Xp is the position of the feature center, in real physical units
    Xp = features.Xcenter .* repmat(features.spacing', ...
				    1,n_features);
    prev_Xp = previous_features.Xcenter .* ...
	      repmat(previous_features.spacing', ...
		     1,size(previous_features.Xcenter,2));
    [sd,closest_nbr_prev] = mindist(Xp,prev_Xp);
    % Assign the offset.
    if params.register_each_stack
      % interpolate in time
      offset_prev = permute(previous_features.registration.offset,[3 1 2]);
      offset_cur = interp1(previous_features.stacknum,...
        offset_prev(:,:,closest_nbr_prev),...
        features.stacknum,...
        'linear',...
        'extrap');
      registration.offset = ipermute(offset_cur,[3 1 2]);
    else
      % Define parameters on just the previously-used stack numbers.
      registration.offset = registration.offset(:,closest_nbr_prev,:);
    end
    % Adjust the scale, to accomodate the changes in pyramid level
    for i = 1:stack_dims
      registration.offset(i,:,:) = registration.offset(i,:,:) * ...
        (previous_features.spacing(i) / features.spacing(i));
    end
    if ~strcmp(registration.model,'translation')
      error('Registration model not implemented');
    end
  else
    % Default behavior when no previous_features available
    registration.model = 'translation';
    registration.offset = zeros(stack_dims,n_features,n_stacks);
    registration.stacknum = features.stacknum;
  end
  features.registration = registration;
  
  
  % Initialize the features themselves, using the image data of the first
  % stack
  Xmin = features.Xmin + repmat(params.margin(:),1,n_features);
  Xmax = features.Xmax - repmat(params.margin(:),1,n_features);
  x = cell(1,stack_dims);
  padx = cell(1,stack_dims);
  for j = 1:stack_dims
    padx{j} = params.margin(j):params.margin(j)+params.center(j)-1;
  end
  pad_image_size = params.center+2*params.margin;
  if stack_dims == 1
    pad_image_size = [1 pad_image_size];
  end
  pad_image = zeros(pad_image_size);
  features.T = nan(n_features,n_stacks);
  features.S = cell(n_features,1);
  for featureIndex = 1:n_features
    for dimIndex = 1:stack_dims
      x{dimIndex} = (Xmin(dimIndex,featureIndex):min(Xmax(dimIndex,featureIndex),features.stacksize(dimIndex)-1))+1;
      padxtmp{dimIndex} = (padx{dimIndex}(1:length(x{dimIndex})))+1;
    end
    pad_image(:) = 0;  % clear any previous data
    pad_image(padxtmp{:}) = im(x{:});
    image_total_intensity = mean(pad_image(:));
    features.T(featureIndex,1:n_stacks) = image_total_intensity; % constant in time
    if (image_total_intensity > 0)
      pad_image = pad_image / image_total_intensity;
    end
    features.S{featureIndex} = pad_image;
  end

function Xo = grid2col(x)
  if (length(x) == 1)
    Xo = x{1}(:);
    return
  end
  X = cell(size(x));
  [X{:}] = ndgrid(x{:});
  for i = 1:numel(X)
    X{i} = X{i}(:);
  end
  Xo = cat(2,X{:});
  