function ssr(smm,stacknums,deformation_filename,svd_filename,options)
% SSR: simultaneous segmentation and registration
% This function is designed to exploit temporal variability in the process of
% registering and segmenting three-dimensional image stacks.  A typical
% application would be imaging neuronal activity, where the change in
% activity over time helps delineate independent structures in the image
% volume.
%
% The approach is to decompose the image into a product of spatial and
% temporal components, and to spatially register the raw data to this
% decomposition.  Mathematically, the approach is to fit a time series in
% the following way:
%     I(g(x,t),t) "=" I_0(x) + sum_i U_i(x) s_i v_i(t)
% where
%    x' = g(x,t) is the deformation of space that brings the volume into
%      alignment at time t,
%    I_0 (optional, see "subtract_mean" below) is the baseline image (the
%      average over time)
%    the sum over i is the sum over components, where the U's are the
%      "image" of the components (normalized), and s_i v_i(t) describes
%      the intensity as a function of time.
%
% This model allows the intensity of individual components to fluctuate
% over time.  This fluctuation, in turn, helps define the spatial extent
% of individual components.  By optimizing the components and the
% deformation g in parallel, the improvement of one leads to an
% improvement in the other.
%
% One important detail: the image volume should be "bricked" in small
% pieces. The idea is to ensure that the number of components needed to
% adequately describe each brick is substantially smaller than the number
% of time points in the time series.  This is necessary if one is to
% derive a representation of the time-series data that is smaller than
% the original data set.
%
% Note that this function, by default, makes a single "pass" through the
% data.  It might be possible to improve further by repeated passes, at
% the cost of additional computation time.
%
% Syntax:
%   ssr(smm,stacknums,deformation_filename,svd_filename,options)
% where
%   smm is a stackmm object (see STACKMM);
%   stacknums is a list containing the stacks that you want to process,
%     in the order in which they are to be processed.
%   deformation_filename is the name of the .mat file to which to save/load
%     information about the deformation applied to each stack (this is
%     updated on each stack to minimize the impact of lost work if a
%     process is terminated)
%   svd_filename is the name of the file in which to save information
%     about the online SVD analysis.
%   options is a structure which, for the first call, must have the
%     following fields set:
%       bricks: a cell array defining the volume bricks (see
%         PREPARE_BRICKS);
%       rmg_options: any options needed for register_multigrid_options to
%         ensure proper registration (include the field "lambda" to set
%         the regularization penalty)
%       n_components: the # of orthogonal components for each brick
%         (e.g., the maximum number of cells expected in a brick volume)
%     The following fields are optional:
%       ratiometric (default true): if true, the first stack is taken as
%         the reference stack, and the PCA is done on the (smoothed, see
%         below) log-ratio of images.
%       smooth_pixels: if supplied (as a vector of # of pixels along each
%         coordinate axis), this will cause the log-ratio to be smoothed before
%         doing PCA.
%       subtract_mean (default false): if true, this performs PCA rather
%         than SVD by subtracting the data mean, thus focusing the SVD on
%         the variance. In ratiometric analysis this is probably of little
%         consequence, because the log-ratio should be centered on 0
%         anyway.
%       register (default true): if false, skips the registration step
%       register_svd (default true): if true, registration is performed to
%         an estimate of the "ideal" data based on a projection onto the
%         SVD. If false, each stack is registered to the base stack.
%       n_vcycles (default 1): the "inner loop" performs an initial
%         registration based on previous results , followed by n_vcycles
%         iterations of (1) projecting the registered image onto the
%         components, (2) computing the "idealized" image from the
%         projected intensities, and (3) registering the data to the
%         idealized image.
%       coords: if supplied, it snips out a region of the stack (i.e.,
%         im(coords{:})).  This can speed up execution considerably.
%       imclass (default 'single'): controls the precision of the
%         registered data.
%       baseline_deformation: if supplied, this acts as a lookup-table
%         for the initial guess for the deformation u (supply a vector of
%         stack #s, containing the stack # to use as the initial
%         guess).
%
% Example:
% This example would process each stack in the entire data set sequentially.
%    smm = stackmm('vno.imagine');
%    szfull = smm.size;
%    % test with a small volume
%    options = struct('coords',{{200:400,300:450,1:40}});
%    imtst = smm(options.coords{:},1);
%    sz = size(imtst);
%    mask = true(sz); mask = trim_mask_edges(mask,[25 25 4]);
%    options.rmg_options = struct('pixel_spacing',[0.7 0.7 5],...
%        'gap_data',3,'lambda',1e4,'mask',mask);
%    bricksz = round(25./options.rmg_options.pixel_spacing); % 25 micron bricks
%    options.bricks = prepare_bricks(sz,bricksz);
%    options.n_components = 10;  % max of 10 cells in a brick
%    ssr(smm,1:szfull(4),'deformation.mat','svd.mat',options)
%
% Alternatively, you could process every 10th stack with a call like
% this:
%    ssr(smm,1:10:sz(4),'deformation.mat','svd.mat',options)
% You could then check the registration quality, and if adequate then go
% back and fill in the rest with a second call:
%    ssr(smm,1:sz(4),'deformation.mat','svd.mat')
% Note that on the second call you do not need to supply options.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 5)
    options = struct;
  end
  options = default(options,'ratiometric',true,'subtract_mean',false,'register',true,'register_svd',true,'n_vcycles',1,'imclass','double');
  first = true;

  n_stacks = length(stacknums);
  base_stack_index = stacknums(1);
  counter = 0;
  err = 0;  % registration error
  if exist(deformation_filename,'file')
    % We are re-starting an existing computation
    deformationdata = load(deformation_filename);
    first = ~isfield(deformationdata,'rmg_params');
    if ~first
      rmg_params = deformationdata.rmg_params;
    end
  else
    % Starting from scratch
    deformationdata = struct('stacknums',[],'u',{{}});
    if options.register
      deformationdata.rmg_options = options.rmg_options;
    end
    if isfield(options,'coords')
      deformationdata.coords = options.coords;
    end
  end
  if exist(svd_filename,'file')
    % We are re-starting an existing computation
    svddata = load(svd_filename);
    % Skip any stacks that have already been processed
    stacknums = setdiff(stacknums,svddata.stacknums);
    base_stack = ssr_loadstack(smm,deformationdata,base_stack_index,options);
  else
    % Starting fresh. Initialize the SVD (use cell arrays because of bricks)
    svddata = struct('stacknums',[],'bricks',{options.bricks},...
      'n_components',options.n_components,'U',{{}},'S',{{}},'V',{{}});
    if options.subtract_mean
      svddata.mean = [];
    end
  end
  
  for stackIndex = stacknums
    % Load the stack data
    stk = ssr_loadstack(smm,deformationdata,stackIndex,options);
    if (stackIndex == base_stack_index)
      base_stack = stk;
      if ~options.register
        base_stack(base_stack == 0) = nan;  % if already registered, encoded as uint16; 0 stands for NaN
      end
    end

    % Find the deformation that should be used as the baseline guess
    u = [];
    if isfield(options,'baseline_deformation')
      % Look up a user-specified value (indexed by absolute stack #)
      baseI = options.baseline_deformation(stackIndex);
      baseII = find(deformationdata.stacknums == baseI);
    else
      % Look up the previously-processed stack that is closest in terms
      % of stack #
      [~,baseII] = min(abs(deformationdata.stacknums-stackIndex));
    end
    if ~isempty(baseII)
      u = deformationdata.u{baseII};
    end

    % Do an initial registration, using our initial guess for the deformation
    if ~isempty(u)
      stkr = register_multigrid_warp(stk,u,rmg_params);
    else
      stkr = stk;
    end
    
    if ~options.register
      stkr(stkr == 0) = nan;  % if already registered, encoded as uint16; 0 stands for NaN
    end
    
    if options.ratiometric
      stkr = ssr_ratio(stkr,base_stack,options);
    end
    
    %
    % Iteratively improve the registration and the predicted image data
    %
    if options.register
      for vcycleIndex = 1:options.n_vcycles
        if options.register_svd && ~isempty(svddata.U)  % only if have components for a prediction...
          % Generate predicted image
          stkr_pred = zeros(size(stkr),options.imclass);
          for brickIndex = 1:length(svddata.bricks)
            % Excise the current brick
            thisbrick = svddata.bricks{brickIndex};
            stkr_snip = stkr(thisbrick{:});
            dataVec = stkr_snip(:);
            if options.subtract_mean
              dataVec = dataVec - svddata.mean{brickIndex};
            end
            thisU = svddata.U{brickIndex};
            thisS = svddata.S{brickIndex};
            % Compensate for NaNs in the registered data
            nanFlag = isnan(dataVec);
            if any(nanFlag)
              % Fill in by smallest deviation from nullspace
              dataVec(nanFlag) = (thisU(nanFlag,:)*thisS) * ((thisU(~nanFlag,:)*thisS)\dataVec(~nanFlag));
            end
            % Project onto the components & calculate predicted brick
            proj = thisU'*dataVec;
            predVec = thisU*proj;
            if options.subtract_mean
              predVec = predVec + svddata.mean{brickIndex};
            end
            predBrick = reshape(predVec,size(stkr_snip));
            % Insert the predicted brick into the stack
            stkr_pred(thisbrick{:}) = predBrick;
          end
          if options.ratiometric
            stkr_pred = exp(stkr_pred) .* base_stack;
          end
        else
          stkr_pred = base_stack;
        end
        if first
          rmg_params = register_multigrid_options(stkr_pred, deformationdata.rmg_options.mask,...
            deformationdata.rmg_options);
          first = false;
          % Compress for disk storage
          deformationdata.rmg_params = register_multigrid_options_compress(rmg_params);
        end
        % Register to this predicted stack
        if (stackIndex ~= base_stack_index)
          rmg_params = register_multigrid_options_expand(rmg_params, stkr_pred);
          [u,err] = register_multigrid_vcycle(u,stk,deformationdata.rmg_options.lambda,rmg_params);
          stkr = register_multigrid_warp(stk,u,rmg_params);
        else
          stkr = stk;
        end
        if options.ratiometric
          stkr = ssr_ratio(stkr,base_stack,options);
        end
      end
    end

    %
    % Update the SVD
    %
    % Prevent regions that are not well-registered from contaminating the
    % SVD
    stkr = stkr .* deformationdata.rmg_options.mask;
    if (~options.ratiometric || stackIndex ~= base_stack_index)
      for brickIndex = 1:length(svddata.bricks)
        % Excise the current brick
        thisbrick = svddata.bricks{brickIndex};
        stkr_snip = stkr(thisbrick{:});
        dataVec = stkr_snip(:);
        inputs = {dataVec,svddata.n_components};  % args for initial call
        if (length(svddata.U) >= brickIndex)
          % This is not the initial call, set up extra input args.
          % If we are re-processing this stack, we need to exclude the old
          % result from the online update
          keepFlag = stackIndex ~= svddata.stacknums;
          svddata.stacknums = svddata.stacknums(keepFlag);
          Vtmp = svddata.V{brickIndex}(keepFlag,:);
          % fixme: if we are replacing, how to fix S?
          inputs(end+1:end+3) = {svddata.U{brickIndex}, ...
            svddata.S{brickIndex},Vtmp};
          if options.subtract_mean
            inputs{end+1} = svddata.mean{brickIndex}; %#ok<AGROW>
          end
          %       else
          %         % This is the initial call. If subtracting the mean, then it makes
          %         % no sense to do SVD on the first call
          %         if (options.subtract_mean)
          %           svddata.mean{brickIndex} = dataVec;
          %           svddata.U{brickIndex} = zeros(length(dataVec),0);
          %           svddata.V{brickIndex} = zeros(0,1);
          %           svddata.S{brickIndex} = zeros(0,0);
          %           continue
          %         end
        end
        if options.subtract_mean
          [svddata.U{brickIndex},svddata.S{brickIndex}, ...
            svddata.V{brickIndex},svddata.mean{brickIndex}] = ...
            crunchy_isvd(inputs{:});
        else
          [svddata.U{brickIndex},svddata.S{brickIndex}, ...
            svddata.V{brickIndex}] = ...
            crunchy_isvd(inputs{:});
        end
      end
      svddata.stacknums(end+1) = stackIndex;
    end
    
    %
    % Save progress to disk
    %
    save(svd_filename,'-struct','svddata');
    if options.register
      % Append/replace deformation data
      storeIndex = find(stackIndex == deformationdata.stacknums);
      if isempty(storeIndex)
        storeIndex = length(deformationdata.stacknums)+1; % append
      end
      deformationdata.stacknums(storeIndex) = stackIndex;
      deformationdata.u{storeIndex} = u;
      deformationdata.err(storeIndex) = err;
      save(deformation_filename,'-struct','deformationdata');
    end
    
    %
    % Show progress to the user
    %
    counter = counter+1;
    fprintf('\nFinished stack %d (%g%% complete)\n\n',stackIndex,100*counter/n_stacks);
  end
end

function stk = ssr_loadstack(smm,deformationdata,stackIndex,options)
  if isfield(deformationdata,'coords')
    stki = smm(deformationdata.coords{:},stackIndex);
  else
    stki = smm(:,:,:,stackIndex);
  end
  stk = cast(stki,options.imclass);
end

function stko = ssr_ratio(stk,base_stack,options)
  rat = stk./base_stack;
  rat(isinf(rat) | isnan(rat)) = 1;
  stko = log(rat);
  if isfield(options,'smooth_pixels')
    nanFlag = isnan(stko);
    stko = imfilter_gaussian(stko,options.smooth_pixels);
    stko(nanFlag) = nan;
  end
end
