function [err,shift] = frame_error(smm,options)
% frame_error: compute mismatch on a per-frame basis for whole movie
%
% Syntax:
%   err = frame_error(smm)
%   [err,shift] = frame_error(smm,options)
% where
%   smm is a stackmm object
%   options is a structure which may have the following fields:
%     register (default true): if true, each stack is first
%       rigid-registered to the "base" stack before computing the framewise
%       mismatch
%     base_stacknum (default mid-movie, just before a stimulus): the stack
%       to use for comparison, when computing the mismatch
%     show_progress (default true): if true, a progress bar is plotted
%     dx_max (default smm.size/2):  limits the fourier method from choosing
%       image edges for translation
%     algorithm (default 'nancorr'):  sets the algorithm for determining 
%        optimal rigid shift to apply to data, note that nancorr requires
%        more fourier transforms (5 to 3) on a nanpadded data set that is 
%        2^3 times larger than phase correlation, resulting in an ~16 fold
%        increase in processing time
% and
%   err is an n_frames-by-n_stacks matrix, giving the mismatch error (mean
%     square pixelwise difference) for each frame in each stack. NaN means
%     the frame was missing due to registration issues
%   shift is an 3-by-n_stacks matrix giving the amount of translation
%     applied to each stack to bring it into register with the base stack.
%
% See also: find_bad_frames.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  header = smm.header;
  sz = smm.size;
  if isempty(header.stim_lookup)
    prestim = 1:sum(header.nstacks);
  else
    prestim = find(header.stim_lookup(1:end-1) == 0 & header.stim_lookup(2:end) > 0);
    prestim = prestim(:);
  end
  options = default(options,'register',true,...
    'base_stacknum',prestim(round(length(prestim)/2)),...
    'show_progress',true,...
    'dx_max',sz(1:3)/2,...
    'algorithm','nancorr');
  if ~options.register && nargout > 1
    error('Can''t request shift if register=false')
  end
  
  fixed = single(smm(:,:,:,options.base_stacknum));
  err = zeros(sz(3:4));
  if options.register
    if strcmp(options.algorithm, 'nancorr')
      % Pre-compute fourier transform of fixed image, to save time
      fixed_padded = register_translate_pad(fixed);
      fixed_fft = register_translate_nancorr(fixed_padded);
    else
      fixed_fft = fftn(fixed);
    end
    shift = zeros(sz(4),3);
    h = fspecial('gaussian',5,1);
    coords = cell(1,3);
  end  
  
%   if options.register
%     shift = zeros(sz(4),3);
%     % Pre-compute fourier transform of fixed image, to save time
%     fixed_fft = fftn(fixed);
%     h = fspecial('gaussian',5,1);
%     coords = cell(1,3);
%   end

  % Prepare for sub-pixel
  rng = cell(1,3);
  minops = optimset('GradObj','on','Display','off');
  if options.show_progress
    pgoptions = struct('max',sz(4));
    tic
  end
  
  for stackIndex = 1:sz(end)
    moving = single(smm(:,:,:,stackIndex));
    if strcmp(options.algorithm, 'nancorr')
      % pad moving image
      moving_padded = register_translate_pad(moving);
      thisshift = register_translate_nancorr(fixed_fft, moving_padded,options);
      imr = image_shift(moving,thisshift);
      shift(stackIndex,:) = thisshift;
    else options.register
      % Register using phase correlation. Could call register_rigid, but
      % that involves an extra fourier transform of the fixed image each
      % time
      moving_fft = fftn(moving);
      rf = moving_fft .* conj(fixed_fft) ./ abs(moving_fft.*fixed_fft);  % FT of normalized phase correlation
      rf(isnan(rf)) = 0;
      r = ifftn(rf);  % normalized phase correlation (for this channel)
      r = fftshift(r);
      r = imfilter(r,h);  % smooth the correlation, so we don't choose something whacky
      rtmp = killedges(r,-inf);  % ensure we don't choose an edge pixel as the max
      [~,indx] = max(rtmp(:));
      [coords{:}] = ind2sub(sz(1:3),indx);
      % Subpixel shift: find the peak using quadratic interpolation
      for dimIndex = 1:3
        rng{dimIndex} = coords{dimIndex} + [-1 0 1];
      end
      rsnip = r(rng{:});
      func = @(g) imqinterp(g,-double(rsnip));
      dx = fmincon(func,[2 2 2],[],[],[],[],1.5*[1 1 1],2.5*[1 1 1],[],minops) - 2; % fractional shift
      thisshift = cat(2,coords{:}) - ceil((sz(1:3)+1)/2) + dx;
      shift(stackIndex,:) = thisshift;
      imr = image_shift(moving,thisshift);      
    end
    mismatch = nanmean(nanmean((fixed-imr).^2,1),2);
    err(:,stackIndex) = mismatch(:);
    if (toc > 3 || stackIndex == sz(end))
      pgoptions.progress = stackIndex;
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end
  if (nargout > 1)
    shift = shift';
  end
  