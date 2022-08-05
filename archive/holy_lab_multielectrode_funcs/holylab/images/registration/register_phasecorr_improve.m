function [uc,imwold] = register_phasecorr_improve(uoldc,fixed,moving,options,imwold)
% register_phasecorr_improve: optimize a deformation for registering images
%
% Syntaxes:
%   u = register_phasecorr_improve(szu,fixed,moving,options)
% This calculates a deformation, with grid size defined by the vector szu,
% given the input images.  "fixed" and "moving" are two (possibly
% multidimensional or multi-color) images.  "moving" means the image you
% want to deform to bring it into alignment with "fixed."
%   As an example, with a two-dimensional image, specifying szu = [1 1]
% would perform a global shift (a rigid deformation).  Specifying
% szu = [3 3] would instead partition the image into 3-by-3 chunks and
% register each chunk separately.  u would be a cell array with two
% entries, each a 3-by-3 array giving the shift (in pixels) applied to
% each chunk.
%   "options" arises from register_phasecorr_intialize.
%
%   [u,imwold] = register_phasecorr_improve(uold,fixed,moving,options)
%   [u,imwold] = register_phasecorr_improve(uold,fixed,moving,options,imwold)
% This syntax improves upon an existing deformation uold.  imwold is the
% moving image warped by the old deformation (uold), which is
% calculated en route to improving the deformation. If you want to
% calculate the moving image deformed by the new solution u, call
% register_phasecorr_warp. Note that if you happen to have already
% calculated the old warped image, you can save some computation time by
% supplying it as the fifth argument.
%
% There is a complete example in the help for register_phasecorr_refine;
% using "fixed", "moving", and "options" as described there, you can use
% this function as follows:
%   u = register_phasecorr_improve([3 3],fixed,moving,options);
%   imw = register_phasecorr_warp(u,moving,options);
%   figure; imshowrgb(fixed,imw)
%
% See also: register_phasecorr_initialize, register_phasecorr_warp,
% register_phasecorr_refine.

% Copyright 2010-2011 by Timothy E. Holy & Julian Meeks

  options = default(options,'lambda',0,'smooth',zeros(1,options.n_dims),'reg_prolong',true);
  
  %% Calculate the "old" warped image
  have_uold = false;
  if iscell(uoldc)
    if isempty(uoldc{1})
      error('If supplying empty u, supply the desired size of u as a vector')
    end
    have_uold = true;
    szu = size(uoldc{1});
    if (nargin < 5)
      imwold = register_phasecorr_warp(uoldc,moving,options);
    end
    uold = cat(options.n_dims+1,uoldc{:});  % Preparation for regularization & composition
  else
    imwold = moving;
    szu = uoldc;
    if (length(szu) ~= options.n_dims)
      error('The size of u does not match the image dimensionality')
    end
    uold = [];
  end
  szu(end+1:options.n_dims) = 1;
   
  %% Compute the phase correlation for image chunks & find the peak in each
  % Generate the grid of cuts. We want the regions to overlap, so define
  % the "left" and "right" edges along each dimension.
  left = cell(1,options.n_dims);
  right = cell(1,options.n_dims);
  edgec = cell(1,options.n_dims);
  for dimIndex = 1:options.n_dims
    n = szu(dimIndex);
    N = options.sz_spatial(dimIndex);
%     edge = (3:2:2*n-1)/(2*(n+1));  % midpoints between evenly-spaced centroids at (1:n)/(n+1), i.e., grid inset inside image
%     f = options.overlap_fraction/(n+1);
    edge = (1:2:2*n-3)/(2*(n-1));  % midpoints between evenly-spaced centroids at (0:n-1)/(n-1), i.e., grid includes corners
    f = options.overlap_fraction/(n-(n>2));
    tmp = ceil([0 edge-f]*(N-1)+1);
    tmp(tmp<1) = 1;
    left{dimIndex} = tmp;
    tmp = floor([edge+f 1]*(N-1)+1);
    tmp(tmp>N) = N;
    right{dimIndex} = tmp;
    edgec{dimIndex} = [1 edge*(N-1)+1 N];  % store for setting limits on search for peak of r
  end
  % Compute the phase correlation for each chunk
  % Initialize storage: for convenience we'll start with the "dimension"
  % being the first coordinate, and later we'll permute so that the
  % dimension coordinate is last.
  u = zeros([options.n_dims szu],options.class);
  coords = cell(1,options.n_dims);
  rng = cell(1,options.n_dims);
  choprng = cell(1,options.n_dims);
  n_chunks = prod(szu);
  if options.lambda > 0
    ucenter = u;
    coef = zeros(szu,options.class);
  end
  for chunkIndex = 1:n_chunks
    [coords{:}] = ind2sub(szu,chunkIndex);
    for dimIndex = 1:options.n_dims
      rng{dimIndex} = left{dimIndex}(coords{dimIndex}):right{dimIndex}(coords{dimIndex});
    end
    r = [];  % this will hold to correlation, averaged across channels
    for valIndex = 1:options.n_values
      fixedsnip = fixed(rng{:},valIndex);
      movingsnip = imwold(rng{:},valIndex);
      if (valIndex == 1)
        % Discard rows, columns, or frames containing NaNs
        % (only need to calculate once, since NaN pattern will be identical
        % across value channels)
        nanrng = nanbox(movingsnip);
        nanrngc = cell(1,size(nanrng,1));  % do this here in case of changed dimensionality
        for dimIndex = 1:size(nanrng,1)
          nanrngc{dimIndex} = nanrng(dimIndex,1):nanrng(dimIndex,2);
        end
      end
      fixedsnip = fixedsnip(nanrngc{:});
      movingsnip = movingsnip(nanrngc{:});
      if isempty(movingsnip)
        break  % will leave r empty
      end
      % edgetaper here?
      f_fft = fftn(fixedsnip);
      m_fft = fftn(movingsnip);
      rf = m_fft .* conj(f_fft) ./ abs(m_fft.*f_fft);  % FT of normalized phase correlation
      rf(isnan(rf)) = 0;
      rchan = ifftn(rf);  % normalized phase correlation (for this channel)
      if isempty(r)
        r = rchan;
      else
        r = r+rchan;
      end
    end
    if isempty(r)
      continue
    else
      r = fftshift(r); % shift so that "zero shift" would be a peak in the center
      ndimsr = ndims(r);  % in rare cases can lose a dimension of r
      if isa(options.smooth,'function_handle')
        r = options.smooth(r);
      elseif any(options.smooth > 0)
        r = imfilter_gaussian_mex(r,options.smooth(1:ndimsr));  % smooth the correlation
      end
      % Find the peak of the cross-correlation
      % Limit the search to points away from regions of overlap
      for dimIndex = 1:options.n_dims
        thisrng = rng{dimIndex}(nanrng(dimIndex,:));  % the range of coordinates of this snippet
        thisedge = edgec{dimIndex}(coords{dimIndex}:coords{dimIndex}+1); % location of edges
        chop = ceil(max(thisedge(1)-thisrng(1),thisrng(2)-thisedge(2))/2); % cut sufficient pixels to stay within edge
        chop = max(chop,0);
        choprng{dimIndex} = chop+1:size(r,dimIndex)-chop;
      end
      rsnip = r(choprng{:});
      if isempty(rsnip)
        continue
      end
      szsnip = size(rsnip);
      ndimsr = ndims(rsnip);
      [~,indx] = max(rsnip(:));
      [coords{1:ndimsr}] = ind2sub(szsnip,indx);
      thisu = cat(2,coords{1:ndimsr}) - ceil((size(rsnip)+1)/2);
      thisu(end+1:options.n_dims) = 0;
      if (options.lambda > 0 && all(thisu == 0))
        % we will do a trial shift anyway to get an estimate of the
        % "stiffness"
        thisu = ones(1,options.n_dims,options.class);
        [coords{:}] = ind2sub(szu,chunkIndex);
        for dimIndex = 1:options.n_dims
          if (coords{dimIndex} > szu(dimIndex)/2)
            thisu(dimIndex) = -1;  % Shift towards center, away from NaNs
          end
        end
      end
      % Calculate how much of a difference this makes, by re-snipping with
      % this shift and calculating the mismatch
      err0 = 0;
      erru = 0;
      rngu = rng;
      for dimIndex = 1:options.n_dims
        rng{dimIndex} = rng{dimIndex}(nanrngc{dimIndex});  % compose the two snipping regions
        tmp = rng{dimIndex} + thisu(dimIndex);  % shift by the desired amount
        nanflag = tmp < 1 | tmp > options.sz_spatial(dimIndex);  % make sure we stay in bounds
        rngu{dimIndex} = tmp(~nanflag);
        rng{dimIndex} = rng{dimIndex}(~nanflag);
      end
      for valIndex = 1:options.n_values
        movingsnip = imwold(rngu{:},valIndex);  % snip with shifted coordinates
        if (valIndex == 1)
          % Clip to avoid NaNs in either the shifted or original moving
          % image
          nanrng = nanbox(movingsnip);
          for dimIndex = 1:size(nanrng,1)
            thisrng = nanrng(dimIndex,1):nanrng(dimIndex,2);
            rngu{dimIndex} = rngu{dimIndex}(thisrng);
            rng{dimIndex} = rng{dimIndex}(thisrng);
          end
          movingsnip = imwold(rngu{:},valIndex);  % re-snip
        end
        fixedsnip = fixed(rng{:},valIndex);
        erru = erru + sum((fixedsnip(:)-movingsnip(:)).^2);
        % Also re-snip unshifted using the same pixels
        movingsnip = imwold(rng{:},valIndex);
        err0 = err0 + sum((fixedsnip(:)-movingsnip(:)).^2);
      end
      if (erru < err0)
          u(:,chunkIndex) = thisu(:);
          if (options.lambda > 0)
            ucenter(:,chunkIndex) = thisu(:);
          end
      end
      if (options.lambda > 0)
        % For this chunk, make the data term be a parabola centered on the
        % best match of the two (comparing shifted vs. unshifted).  Coef is
        % the coefficient of the quadratic, set so that it gives the value
        % err0 at u=0 and erru at the shifted u.
        coef(chunkIndex) = abs(erru-err0)/sum(thisu.^2);
      end
    end
  end
  % Now make the dimension coordinate last
  permorder = [2:options.n_dims+1 1];
  u = permute(u,permorder);
  if options.lambda > 0
    ucenter = permute(ucenter,permorder);
  end
  chunksz = reshape(options.sz_spatial./szu,[ones(1,options.n_dims) options.n_dims]);
  
  %% If necessary, optimize the placements respecting regularization
  if (options.lambda > 0 && have_uold && ~all(szu == 1))
    % Convert the data-term penalty to a mean over pixels, so lambda
    % has more meaningful units
    coef = coef/prod(options.sz_spatial);
    % Compute the part of the penalty that arises from composition with
    % uold
    detJprev = register_logdetpenalty(bsxfun(@rdivide,uold,chunksz));
    if any(detJprev(:) < 0)
      error('The previous deformation had a tear (detJ < 0)');
    end
    % Prolong detJprev, because the penalty will be evaluated on the
    % next-finer grid (see rpi_regularization for discussion). We have to
    % pad it by one entry in all dimensions to make prolongation work, and
    % then snip one off.
    if options.reg_prolong
      szuh = register_phasecorr_prolong(uold,options);
      if isscalar(detJprev)
        dJp = ones(szu,options.class);
      else
        dJp = detJprev;
        colons = repmat({':'},1,options.n_dims);
        for dimIndex = 1:options.n_dims
          thiscol = colons;
          thiscol{dimIndex} = szu(dimIndex);
          dJp(thiscol{:}) = 1;  % 1 is the expected value
        end
      end
      dJph = array_prolong(dJp,szuh);
      for dimIndex = 1:options.n_dims
        rng{dimIndex} = 1:szuh(dimIndex)-1;
      end
      dJph = dJph(rng{:});
    else
      dJph = detJprev;
    end

    % Ensure that the regularization penalty is not infinite
    penalty = @(u) rpi_penalty(u,ucenter,coef,options,chunksz,dJph);
    for iter = 1:10
      [~,val] = rpi_regularization(u,chunksz,dJph,options);
      if (~isinf(val))
        break
      end
      % Average u with a smoothed version of u
      warning('register:phasecorr','Smoothing u');
      u = (u + array_prolong(array_restrict(u),szu))/2;
    end
    if (isinf(val))
      error('Could not make regularization penalty finite');
    end
    
    u = conjgrad(penalty,u,struct('iter_max',5000));
  end
  
  %% Compose with the previous solution to get the improved deformation
  uc = register_phasecorr_composeu(uold,u,chunksz);
end

%% Penalty function for regularization
function [val,grad] = rpi_penalty(u,ucenter,coef,options,chunksz,detJprev)
  % Compute the data portion
  du = u-ucenter;
  val = sum(du.^2,options.n_dims+1).*coef;
  val = sum(val(:));
  grad = bsxfun(@times,2*coef,du);
  % Add the regularization portion
  szu = size(u);
  if (options.lambda > 0 && any(szu(1:end-1) > 1))
    [~,val_reg,grad_reg] = rpi_regularization(u,chunksz,detJprev,options);
    val = val + options.lambda*val_reg;
    grad = grad + options.lambda*grad_reg;
  end
end

function [detJnew,val_reg,grad_reg] = rpi_regularization(u,chunksz,detJprev,options)
  du = bsxfun(@rdivide,u,chunksz);
  % Evaluate penalty on a finer grid, so future u prolongation doesn't
  % cause trouble. (Prolonging u and calculating detJ does not always
  % yield results that are terribly close to calcuating detJ and then
  % prolonging detJ.)
  if options.reg_prolong
    [~,duhc] = register_phasecorr_prolong(du,options);
    duh = cat(options.n_dims+1,duhc{:});
  else
    duh = du;
  end
  if isempty(detJprev)
    [detJnew,val_reg,grad_reg] = register_logdetpenalty(duh);
  else
    % detJprev was already prolonged to a finer grid
    [detJnew,val_reg,grad_reg] = register_logdetpenalty(duh,1,detJprev);
  end
  if options.reg_prolong
    rflag = size(duh) > size(du);
    grad_reg = 2^sum(rflag)*array_restrict(grad_reg,rflag);
  end
end