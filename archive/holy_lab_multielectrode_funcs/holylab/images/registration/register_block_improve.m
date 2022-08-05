function [uc,imwold,options] = register_block_improve(ucold,fixed,moving,options,imwold)
  % register_block_improve: calculate or refine a deformation to improve image
  % registration
  %
  % syntax:
  %   uc = register_block_improve(usz,fixed,moving,options)
  % This syntax calculates a deformation to align a reference image
  % ("fixed") and a warped image ("moving").
  % where
  %   usz is a vector specifying the dimensions of the grid defining the
  % deformation. options is set by call register_phasecorr_initialize, with some
  % additional settings. See the example below.
  %
  %   uc = register_block_improve(ucold,fixed,moving,options)
  % This syntax improves upon an original deformation ucold.
  %
  % Example:
  % % Create fixed and moving images
  % im = double(imread('cameraman.tif'));
  % scale = 10;  % max # of pixels in shift
  % sz = size(im);
  % udef = randn([sz 2]);
  % udef = imfilter_gaussian(udef,[30 30 0]);  % smooth the random deformation
  % udef = udef * (scale / max(abs(udef(:)))); % set the scale
  % imw_full = register_multigrid_warp(im,udef);
  % rng = 15:240;
  % moving = imw_full(rng,rng);
  % fixed = im(rng,rng);
  % figure; imshowrgb(fixed,moving); title('Unregistered images');
  %
  % % Define the registration parameters
  % pyramid = array_restrict_schedule(size(fixed));
  % ops = register_phasecorr_initialize(fixed,struct('pyramid',pyramid));
  % ops.lambda = 1e4;
  % ops.plot = true;
  %
  % % Do the registration
  % % (note [5 5] is one of the sizes in the pyramid)
  % usz = [5 5];
  % uc = register_block_improve(usz,fixed,moving,ops);
  %
  % % Show the results
  % imw = register_phasecorr_warp(uc,moving,ops);
  % figure; imshowrgb(fixed,imw); title(sprintf('[%d %d]',size(uc{1})));
  %
  % % Also show the deformation itself
  % funch = register_gui_utilities;
  % unorm = funch.upc2mg(uc,ops.sz_spatial);
  % figure; imflow_display(unorm);
  
  % Copyright 2011 by Timothy E. Holy
  
  
  colons = repmat({':'},1,options.n_dims);
  
  %% Calculate the "old" warped image
  if iscell(ucold)
    if isempty(ucold{1})
      error('If supplying empty u, supply the desired size of u as a vector')
    end
    usz = size(ucold{1});
    if (nargin < 5)
      imwold = register_phasecorr_warp(ucold,moving,options);
    end
  else
    imwold = moving;
    usz = ucold;
    if (length(usz) ~= options.n_dims)
      error('The size of ucold does not match the image dimensionality')
    end
    ucold = [];
  end
  usz((end+1):options.n_dims) = 1;
  
  %% Set defaults
  initial_guess_filter = repmat(3,1,options.n_dims);
  initial_guess_filter(usz < initial_guess_filter) = 1;
  options = default(options,...
    'descent','greedy',...
    'initial_guess','global',...
    'initial_guess_filter',initial_guess_filter,...
    'smooth',[],...
    'lambda',0,...
    'barrier_mu',-1.0);
  
  %% Compute the mismatch
  if length(usz) == 2
    mismatch = register_block_mismatch(usz,fixed,imwold,options);
  elseif length(usz) == 3 % use cuda accelerated code
    mismatch = register_block_mismatch_cuda(usz,fixed,imwold,options);
  end
  n_blocks = numel(mismatch);
  % If requested, smooth the mismatch
  if ~isempty(options.smooth)
    if isa(options.smooth,'function_handle')
      for blockIndex = 1:n_blocks
        mismatch{blockIndex} = options.smooth(mismatch{blockIndex});
      end
    elseif any(options.smooth > 0)
      for blockIndex = 1:n_blocks
        mismatch{blockIndex} = imfilter_gaussian_mex(mismatch{blockIndex},...
          options.smooth(1:ndims(mismatch{blockIndex})));
      end
    end
  end
  
  % Set the edges to max(mismatch{blockIndex}) so they will be avoided
  for blockIndex = 1:n_blocks
    mismatch{blockIndex} = killedges(mismatch{blockIndex},...
      max(mismatch{blockIndex}(:)));
    % Padding the edge to Inf sometimes will cause a error reporting
    % requiring a not Inf initial value, especially for 2D image
    % registration.
%     mismatch{blockIndex} = killedges(mismatch{blockIndex},...
%       Inf);
  end
  
  %% Set up the initial guess
  centerc = cell(1,options.n_dims);
  uc = cell(1,options.n_dims);
  coords = cell(1,options.n_dims);
  for dimIndex = 1:options.n_dims
    centerc{dimIndex} = zeros(usz,options.class);
    uc{dimIndex} = zeros(usz,options.class);
  end
  for blockIndex = 1:n_blocks
    this_sz = size(mismatch{blockIndex});
    this_center = ceil(this_sz/2);
    this_center((end+1):options.n_dims) = 1;
    for dimIndex = 1:options.n_dims
      centerc{dimIndex}(blockIndex) = this_center(dimIndex);
    end
    switch options.initial_guess
      case 'current_position'
        % Do nothing, because the correct initial offset in uc is zero
      case 'global'
        % Set the initial offset to be at the minimum of the mismatch for
        % this block
        [~,minIndex] = min(mismatch{blockIndex}(:));
        [coords{1:length(this_sz)}] = ind2sub(this_sz,minIndex);
        for dimIndex = 1:length(this_sz)
          uc{dimIndex}(blockIndex) = coords{dimIndex} - this_center(dimIndex);
        end
      otherwise
        error('initial_guess setting not recognized');
    end
  end
  
  %% debug code
%   minSum = 0.0;
%   for blockIndex = 1:n_blocks
%     minSum = minSum + min(mismatch{blockIndex}(:));
%   end
  
  if ~isempty(options.initial_guess_filter)
    for dimIndex = 1:options.n_dims
      if ndims(uc{dimIndex}) == 2
        uc{dimIndex} = medfilt2(uc{dimIndex},options.initial_guess_filter);
      elseif ndims(uc{dimIndex}) == 3
        uc{dimIndex} = medfilt3(uc{dimIndex},options.initial_guess_filter);
      end
    end
  end
  
  %% Establish the geometry of the u-grid
  %   corner = cell(1,options.n_dims);
  % %   interior = cell(1,options.n_dims);
  %   for dimIndex = 1:options.n_dims
  %     n = usz(dimIndex);
  %     corner{dimIndex} = 1:n;
  % %     interior{dimIndex} = (n-1)/(n+1)*(1:n)+1;
  %   end
  
  %% Prepare the interpolation coefficients
  % The composition will be calculated by quadratic interpolation; to
  % ensure accurate values, we need to first compute interpolation
  % coefficients.
  ucoldcoef = ucold;
  if ~isempty(ucold)
    for dimIndex = 1:options.n_dims
      ucoldcoef{dimIndex} = qinterp_grid_inverse(ucold{dimIndex},'reflect');
    end
  end
  
  %% Set the penalty function
%   func = @(u) register_block_penalty(mismatch,ucold,u,interior,options);
%   func = @(u) register_block_penalty(mismatch,ucoldcoef,u,options);
  func = @(u,barrier_mu) register_block_penalty(...
    mismatch,ucoldcoef,u,barrier_mu,options);
  % FIXME: do something to ensure that the initial regularization penalty is not
  % infinite
  
  %% Minimize the penalty
  u = cat(options.n_dims+1,uc{:});
  switch options.descent
    case 'greedy'
      ops = options;
      ops.Display = true;
      ops.iter_max = 10*numel(u);
      [u,fval] = conjgrad(func,u,ops);
      % determine whether register computation is sufficient or not
%       if options.fval > fval(1) 
%         options.fval = fval(1);
%       else
%         options.isdone = true;
%       end
      
      % u = fminunc(func,u(:),optimset('Display','iter','GradObj','on','DerivativeCheck','off'));
      % u = reshape(u,[usz options.n_dims]);
    otherwise
      error('descent setting not recognized');
  end
  for dimIndex = 1:options.n_dims
    uc{dimIndex} = u(colons{:},dimIndex);
  end
  
  %% Compose old and new
  %   uc = rbi_compose(ucold,uc,interior,options);
  %   if ~isempty(ucold)
  %     ucold{:}
  %     ucoldcoef{:}
  %   end
  %   uc{:}
  uc = rbi_compose(ucoldcoef,uc,options);
 
  % Plotting. Delete from production code ->
%   funch = register_gui_utilities;
%   unorm = funch.upc2mg(uc,options.sz_spatial);
%   figure; imflow_display(unorm); title('Initial u');
  % <-   
  
function ucComp = rbi_compose(ucOld,ucNew,options)
  % This uses register_logdetpenalty_composition to perform the
  % composition, so that there is no risk of disagreement
  gridsz = size(ucNew{1});
  blocksz = options.sz_spatial./gridsz;
  % Scale the deformation to a unit-pixel grid
  if ~isempty(ucOld)
    for dimIndex = 1:options.n_dims
      ucOld{dimIndex} = ucOld{dimIndex}/blocksz(dimIndex);
    end
  end
  for dimIndex = 1:options.n_dims
    ucNew{dimIndex} = ucNew{dimIndex}/blocksz(dimIndex);
  end
  ucComp = register_logdetpenalty_composition(ucOld,ucNew);
  % Undo the scaling
  for dimIndex = 1:options.n_dims
    ucComp{dimIndex} = ucComp{dimIndex}*blocksz(dimIndex);
  end
