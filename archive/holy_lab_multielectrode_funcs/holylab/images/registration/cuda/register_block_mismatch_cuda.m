function mismatch = register_block_mismatch_cuda(usz,fixed,moving,options)
% register_block_mismatch_cuda: compute image mismatch in tapered blocks
%
% This function is the same as register_block_mismatch() with the only
% difference that this function will load the cuda acclerated
% register_mismatch_nancorr_cuda().
%
% This routine splits fixed and moving images into blocks, and computes the
% mismatch (the sum of square differences) in each block. Each block is
% tapered with a separable triangular window that, when summed over blocks,
% sums to 1 at each pixel. Using Fourier methods, the mismatch is computed
% for all integer-pixel translations of the moving image block (up to a
% maximum shift, which defaults to half the block size).
%
% The main sophistication behind the splitting lies in the marking of
% pixels for which data are not available---handling "missing data." See
% register_mismatch_nancorr. Since Fourier methods are used, all blocks are
% zero-padded to ensure that translations are not contaminated by periodic
% boundary conditions.
%
% Syntax:
%   mismatch = register_block_mismatch_cuda(usz,fixed,moving,options)
% where
%   usz is the number of blocks along each coordinate
%   fixed is the fixed image
%   moving is the moving image
%   options comes from register_phasecorr_initialize. You can add extra
%     fields:
%       maxshift: a 1-by-n_dims vector (or optionally a scalar), specifying
%         the maximum amount of translation (in pixels) allowed between the
%         blocks of the fixed and moving images.
%       plot: if true, information about each block will be displayed.
% and
%   mismatch is a cell array of size usz, where each element is a numeric
%     array returning the mismatch as a function of displacement.  The
%     center of each array corresponds to zero displacement.
%
% See also: register_mismatch_nancorr.

% Copyright 2011 by Timothy E. Holy

% Fixme: implement RGB/multi-valued images

  %% Argument parsing
  options = default(options,'plot',false);
  have_maxshift = isfield(options,'maxshift');
  N = prod(options.sz_spatial);
  %% Define the block boundaries
  center = cell(1,options.n_dims);
  left = cell(1,options.n_dims);
  right = cell(1,options.n_dims);
  for dimIndex = 1:options.n_dims
    if usz(dimIndex) > 1
      tmp = round(linspace(1,options.sz_spatial(dimIndex),usz(dimIndex)));
      center{dimIndex} = tmp;
      left{dimIndex} = tmp([1 1:end-1]);
      right{dimIndex} = tmp([2:end end]);
    else
      center{dimIndex} = ceil(options.sz_spatial(dimIndex)/2);
      left{dimIndex} = 1;
      right{dimIndex} = options.sz_spatial(dimIndex);
    end
  end
  %% Allocate storage and initialize any variables
  mismatch = cell(usz);
  n_blocks = prod(usz);
  x_src = cell(1,options.n_dims);
  x_dest = cell(1,options.n_dims);
  span_from_fixed = zeros(options.n_dims,2);
  n = zeros(1,options.n_dims);
  [x_src{:}] = ind2sub(usz,(1:n_blocks)');  % just using x_src temporarily
  coord_all = cat(2,x_src{:});
  %% Process each block
  for blockIndex = 1:n_blocks
    coord = coord_all(blockIndex,:);
    %% Snip out block from the fixed image, and build the taper (weight) w
    w = 1;  % tapering function
    for dimIndex = 1:options.n_dims
      thiscoord = coord(dimIndex);
      thisspan = [left{dimIndex}(thiscoord) center{dimIndex}(thiscoord) right{dimIndex}(thiscoord)];% coordinates to be snipped
      span_from_fixed(dimIndex,:) = thisspan([1 3]);
      x_src{dimIndex} = thisspan(1):thisspan(3);
      n(dimIndex) = length(x_src{dimIndex});
      if (n(dimIndex) == 2)
        error('Cannot make blocks smaller than 3'); % except singleton dimensions
      end
      if (n(dimIndex) == 1)
        w1 = 1;
      else
        dn = 1./diff(thisspan);
        w1 = [0:dn(1):1-dn(1) 1 1-dn(2):-dn(2):0];
        w1sz = ones(1,options.n_dims);
        w1sz(dimIndex) = n(dimIndex);
        w1 = reshape(w1,w1sz);
      end
      w = bsxfun(@times,w,w1);
    end
    fixedsnip = fixed(x_src{:});
    if options.plot
      subplot(1,2,1); imshowrgb(fixedsnip,moving(x_src{:})); 
      title(sprintf('block %d',blockIndex));
    end
    %% Determine the size of the zero-padded array
    if have_maxshift
      maxshift = (n-1)/options.maxshift;
    else
      maxshift = (n-1)/6;
    end
    maxshift = floor(maxshift);
    fftsz = 2*n + 2*maxshift; % size of the FFT transform
    % round up on a 32 base, speeding up of FFTN transform
    fftsz = ceil(fftsz/32) * 32;
    
    %% Determine the snip coordinates of the moving image
    % Note we will snip out a larger region of the moving image so that
    % we have actual pixel data for all the shifts
    span_of_moving = bsxfun(@plus,span_from_fixed,maxshift(:)*[-1 1]);
    span_from_moving = span_of_moving;
    span_from_moving(span_from_moving<1) = 1;
    flag = span_from_moving(:,2) > options.sz_spatial(:);
    span_from_moving(flag,2) = options.sz_spatial(flag)';
    %% Embed the fixed snippet and the weight
    fixedpad = zeros(fftsz,options.class);
    wpad = zeros(fftsz,options.class);
    ds = span_from_fixed(:,1) - span_of_moving(:,1);
    for dimIndex = 1:options.n_dims
      thisspan = ([1 n(dimIndex)])+ds(dimIndex);
      x_dest{dimIndex} = thisspan(1):thisspan(2);
    end
    fixedpad(x_dest{:}) = fixedsnip;
    wpad(x_dest{:}) = w;
    %% Embed the moving block in a zero-padded array
    movingpad = zeros(fftsz,options.class);
    % Missing data will be marked with NaN, let's pre-mark it now
    ds = diff(span_of_moving,1,2)+1;
    for dimIndex = 1:options.n_dims
      thisspan = [1 ds(dimIndex)];
      x_dest{dimIndex} = thisspan(1):thisspan(2);
    end
    movingpad(x_dest{:}) = NaN;
    % Embed the available data from the moving image
    ds = span_from_moving(:,1) - span_of_moving(:,1);
    for dimIndex = 1:options.n_dims
      thisspan = [1 diff(span_from_moving(dimIndex,:))+1]+ds(dimIndex);
      x_dest{dimIndex} = thisspan(1):thisspan(2);
      x_src{dimIndex} = span_from_moving(dimIndex,1):span_from_moving(dimIndex,2);
    end
    movingpad(x_dest{:}) = moving(x_src{:});
    
    %% check the padded image block is not too big
    if((numel(fixedpad)*8/1024^3) > 0.7)
      error('The padded image block is too big. Use BIGGER USZ instead.');
    end
    
    %% debug code, calculate the total memory for saving FFT
%     outfile = 'output.out';
%     fid = fopen(outfile,'r');
%     previous = fscanf(fid,'%g',[1]);
%     if(isempty(previous))
%       previous = 0.0;
%     end
%     fclose(fid);
%     fid = fopen(outfile,'w');
%     current = previous + numel(fixedpad)*8/1024^3*3;
%     fprintf(fid, '%g', current);
%     fprintf(fid, '\n');
%     fclose(fid);
    
    
    %% Compute the correlation
    [numerator,denominator] = ...
        register_mismatch_nancorr_cuda(fixedpad,movingpad,wpad);
    
    %% Snip out a region of size maxshift on either side of the center
    thiscenter = ceil((fftsz+1)/2);
    for dimIndex = 1:options.n_dims
      thisspan = thiscenter(dimIndex) + [-1 1]*maxshift(dimIndex);
      x_src{dimIndex} = thisspan(1):thisspan(2);
    end
    numerator = numerator(x_src{:});
    if ~isscalar(denominator)
      denominator = denominator(x_src{:});
    end
    mismatch{blockIndex} = numerator./denominator/N;
    %eps_denominator = sqrt(eps(class(denominator)))*max(denominator(:));
    %mismatch{blockIndex}(denominator < eps_denominator) = inf;
    
    if options.plot
      [~,minIndex] = min(mismatch{blockIndex}(:));
      msz = size(mismatch{blockIndex});
      [x_src{:}] = ind2sub(msz,minIndex);
      coord = [x_src{:}] - ceil(msz/2);
      subplot(1,2,2); imagesc(mismatch{blockIndex}); axis image; 
      title(sprintf('%d ',coord)); pause
    end
    
  end
   
