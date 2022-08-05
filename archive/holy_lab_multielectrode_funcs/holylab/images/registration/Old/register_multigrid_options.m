function options = register_multigrid_options(imFixed,options)
  if (nargin < 2)
    options = struct;
  end

  im_sz = size(imFixed);
  im_dimIndex = im_sz > 1;
  im_sz = im_sz(im_dimIndex); % eliminate unity dimensions
  n_dims = length(im_sz);
  n_dims_orig = length(im_dimIndex);
  
  options.n_dims = n_dims;
  options = default(options,'covariant',true);
  default_pixel_spacing = zeros(1,n_dims_orig); default_pixel_spacing(im_dimIndex) = 1;
  options = default(options,'pixel_spacing',default_pixel_spacing);
  options = default(options,'min_pixels',3);
  options = default(options,'relaxation_iterations',[1 1]); % [down up] V-cycle
  options = default(options,'first_relaxation_iterations',options.relaxation_iterations);
  options = default(options,'max_iter',20);  % total # of V-cycle iterations
  options = default(options,'max_coarsest_iter',10); % for solve on coarsest
  
  % Generate the image grid pyramid
  % The main idea is to decimate in a way that brings the pixel spacing
  % into approximate correspondence: if pixels are spaced more closely
  % along some dimensions than others, then decimate along those
  % dimensions first until the spacing is similar in all directions.
  pixel_spacing = options.pixel_spacing(im_dimIndex);
  sz = im_sz;
  image_grid = struct('sz',size(imFixed),...
		      'decimate',zeros(1,n_dims_orig),...
		      'pixel_spacing',options.pixel_spacing);
  min_spacing = min(pixel_spacing);
  decimation_flag = pixel_spacing < sqrt(2)*min_spacing; % dims to decimate
  %sz_break_flag = sz > options.min_pixels & sz./2.^decimation_flag < options.min_pixels;
  sz_break_flag = (sz-decimation_flag)./2.^decimation_flag < options.min_pixels;
  %while all((sz-1)/2 >= options.min_pixels)
  while any(decimation_flag & ~sz_break_flag) && ~any(decimation_flag & sz_break_flag)
    sz_fac = 2.^decimation_flag;
    sz = floor((sz-decimation_flag)./sz_fac);  % toss the odd pixel at end
    pixel_spacing = pixel_spacing .* sz_fac;
    image_grid(end+1).sz = ones(1,n_dims_orig);
    image_grid(end).sz(im_dimIndex) = sz;
    image_grid(end).decimate = false(1,n_dims_orig);
    image_grid(end).decimate(im_dimIndex) = decimation_flag;
    image_grid(end).pixel_spacing = zeros(1,n_dims_orig);
    image_grid(end).pixel_spacing(im_dimIndex) = pixel_spacing;
    % The last items are there to add some fields ahead of time, so that
    % later stages don't create incompatible structures by adding new fields
    image_grid(end).sqrtdetJ0 = [];
    image_grid(end).mu = [];
    min_spacing = min(pixel_spacing);
    decimation_flag = pixel_spacing < sqrt(2)*min_spacing; % dims to decimate
    %sz_break_flag = sz > options.min_pixels & sz./2.^decimation_flag < options.min_pixels;
    sz_break_flag = (sz-decimation_flag)./2.^decimation_flag < options.min_pixels;
  end
  n_grids = length(image_grid);
  
  % Pre-calculate the coarser versions of imFixed
  if options.covariant
    image_grid(1).psiF = sqrt(imFixed);
  else
    image_grid(1).psiF = imFixed;
  end
  for gridIndex = 2:n_grids
    imFixed = imrestrict(imFixed,image_grid(gridIndex).decimate);
    if options.covariant
      image_grid(gridIndex).psiF = sqrt(imFixed);
    else
      image_grid(gridIndex).psiF = imFixed;
    end
  end
  
  options.image_grid = image_grid;
    