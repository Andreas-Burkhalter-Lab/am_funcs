function options = register_multigrid_options(imFixed,mask,options)
% REGISTER_MULTIGRID_OPTIONS: initialize multigrid registration
%
% This function is called before you do any multigrid registration, to
% set the parameters that will control the algorithm. Among other
% settings, it generates a sequence of coarser-resolution versions of the
% fixed image.  If you are registering multiple images (or stacks) to a
% common fixed image, you only need to call this once.
%
% Syntax:
%   options = register_multigrid_options(imFixed,mask)
%   options = register_multigrid_options(imFixed,mask,default_options)
% where
%   imFixed is the fixed image (can be arbitrary-dimensional; if one
%     dimensional, make sure it is a column vector)
%   mask is a logical array of the same size as imFixed; it should be
%     true for the pixels you want to enforce registration for
%     (typically, set to false around the edges, see TRIM_MASK_EDGES)
%   default_options are fields that you want to specify manually (see
%     below for a description).
%
% The output (and the default_options input) is a structure, for which
% the most important fields are the following:
%   pixel_spacing: a vector giving the spacing between adjacent pixels
%     along each coordinate.  For example, if your pixels are spaced 0.25
%     microns in x and y, and 5 microns in z, this would be [0.25 0.25
%     5]. Default: all 1s.
%     Note #1: there is no consequence to changing this in the output
%     of this function; it affects the results only if you specify it as
%     part of the default_options input.  For images with three
%     dimensions it is frequently quite important that you specify this
%     field.
%     Note #2: the process of making coarser-resolution images will
%     selectively restrict along dimensions that are of substantially-higher
%     (>sqrt(2)-fold) resolution, so that as you progress to coarser grids
%     all of the dimensions ultimately have similar pixel spacing.  If you
%     want to bias this process (for example, if your volume is much
%     smaller, in physical space, along one axis than any of the others),
%     you can safely "lie" about this so that your coarser-resolution arrays
%     do not "collapse" along the short dimension.
%   gap_data: an integer that controls the number of degrees of freedom
%     in the deformation, conventionally denoted "u".  gap_data = 0
%     corresponds to "full resolution", where each pixel has its own
%     independent "lookup coordinate." You probably don't want to do
%     this.  Instead, you probably want to define u on a coarser grid and
%     prolong it up to full size, so that u has fewer degrees of freedom
%     than there are pixels in the image.  The default setting is chosen
%     to permit as many grids as possible; to coarsen u further, increase
%     gap_data.
%   wcycle (default true): if true, register_multigrid_vcycle will
%     perform a w-cycle rather than a v-cycle.  A w-cycle takes longer,
%     but may give you a better match.  See any good text on multigrid
%     methods for an explanation of v-cycles and w-cycles; briefly, a
%     v-cycle goes down to the coarsest resolution and directly back up
%     to fine, whereas a w-cycle "staggers" back up, performing a "short"
%     v-cycle after climbing each level of the grid.
%   options.zero_nans (default true): if true, pixels that are NaN after
%     interpolation are set to zero.  Setting this to true can allow the
%     algorithm to "limp along" better, at the potential cost of not
%     diagnosing a problem that could be fixed by making the mask smaller.
%
% Other fields of lesser importance are:
%   n_grids (not present by default): if you find that the algorithm is
%     not making much progress at the lowest grid levels, you can set
%     this field to prevent it from descending too far.  The easiest
%     approach is to examine the command-window output telling you how
%     many grids it would plan to use, and then choosing a lower number.
%     However, before doing this you should strongly consider just
%     increasing gap_data (see above).  Perhaps one exception to this
%     recommendation is if you are polishing a deformation by using
%     multiple v/w-cycles.
%   n_relaxations (default 1): controls the number of "relaxation"
%     operations performed after rising each level of the grid
%     hierarchy.  Increasing this makes the algorithm slower but may
%     decrease your error; weigh this against the benefits of performing
%     another v/w-cycle.
%   min_pixels (default 7): images with fewer than this number of "good"
%     (non-NaN) pixels along each axis are not used.  With
%     this field and the next: if you are unhappy with the way the
%     image_grid (see below) turned out, your best option is to play with
%     pixel_spacing instead.
%   min_g_pixels (default 2): controls the smallest number of pixels
%     along each axis of u; there will be no coarser grids once one of
%     the axes gets to this point.
%   detJvalue (default 1): set this to something other than 1 if you want
%     to stretch or shrink the moving image to fit the fixed image
%   display (default 'all'): set to 'all' (most output in the command
%     window), 'levels' (intermediate), or 'none' (no output).
% There are yet other fields, but these are largely relics.
%  
% An important thing to know: especially if you are a new user or are
% tackling a new data set, you should examine the "image_grid" field
% (itself a structure array) of the output.  There is a great deal of
% information about how it goes from fine-to-coarse (or vice versa) in this
% field.  You can see the fields by doing this:
%    rmg_params.image_grid
% and, for example, get the sizes of all the coarser grids by doing this:
%    sz = {rmg_params.image_grid.sz};
% While this defines the grid hierarchy, you may find that the imFixed
% field of the coarsest grids is empty; this means that the image did not
% have at least min_pixels valid pixels along at least one of the axes.
% See comments above.  These "empty" elements of the grid hierarchy
% nevertheless contain important information about coarse grids for u,
% which typically goes to coarser grids than does imFixed.
%
% See also: REGISTER_MULTIGRID_VCYCLE, REGISTER_MULTIGRID_OPTIONS_COMPRESS,
% REGISTER_DEMO, ARRAY_RESTRICT_SCHEDULE.
  
% Copyright 2009-2010 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end

  im_sz = size(imFixed);
  dimFlag = im_sz > 1;
  if ~dimFlag(1)
    error('First dimension must not be unity (try supplying the transpose image)');
  end
  im_sz = im_sz(dimFlag); % eliminate unity dimensions
  n_dims = length(im_sz);
  n_dims_orig = length(dimFlag);
  imclass = class(imFixed);
  
  options.n_dims = n_dims;
  default_pixel_spacing = zeros(1,n_dims_orig); default_pixel_spacing(dimFlag) = 1;
  options = default(options,'pixel_spacing',default_pixel_spacing);
  options = default(options,'options',struct);
  options.options = default(options.options,'zero_nans',true);
  options = default(options,'min_pixels',7);  % need 7 non-NaN pixels (buffer of 2 per edge -> 3 "good" pixels)
  options = default(options,'min_g_pixels',2);  % relevant if g_level_gap > 0
  options = default(options,'n_relaxations',1); % # of gradient-descent steps after coarse-grid correction
  options = default(options,'max_vcycles',20);  % total # of V-cycle iterations
  options = default(options,'detJvalue',1);  % the desired scaling factor
  if ~isfield(options,'gap_data')
    gapmode = 'auto';
  else
    gapmode = 'manual';
  end
  options = default(options,'gap_data_mode',gapmode); % choices: 'auto','manual'
  options = default(options,'gap_data',2); % the # of restrictions/prolongations separating g and the image data
  options = default(options,'gap_regularize',0); % determines whether g regularization occurs on full grid or native (smaller) grid (only relevant if gap > 0)
  options = default(options,'display','all'); % choices: 'all','levels','none'
  options = default(options,'max_coarsest_iter_mode','manual');  % choices: 'auto','manual'
  options = default(options,'max_coarsest_iter',100); % for conjgrad on coarsest grid
  options = default(options,'cg_fval_tol',1e-4); % tolerance on conjgrad line minimizations
  options = default(options,'linemin_fval_tol',0.01);
  options = default(options,'coarsest_fval_tol',0.01);
  options = default(options,'wcycle',true);
  %options = default(options,'TolX',1e-2);  % tolerance on line minimzations
  mask = logical(mask);

  % Generate the image grid pyramid
  % The main idea is to decimate in a way that brings the pixel spacing
  % into approximate correspondence: if pixels are spaced more closely
  % along some dimensions than others, then decimate along those
  % dimensions first until the spacing is similar in all directions.
  arsoptions = struct('min_pixels',min(options.min_pixels,options.min_g_pixels),...
    'pixel_spacing',options.pixel_spacing,'truncate',true);
  image_grid = array_restrict_schedule(size(imFixed),arsoptions);
  
  % Fill in the image grid with data
  g0c = register_g0(im_sz,imclass);
  g0 = cat(n_dims+1,g0c{:});
  imFixed = imqinterp(g0,imFixed);  % provide the same smoothing imMoving will get
  if options.options.zero_nans
    nanFlag = isnan(imFixed);
    imFixed(nanFlag) = 0;
  end
  image_grid(1).imFixed = imFixed;
  image_grid(1).mask = mask;
  image_grid(1).g0 = g0;
  
  for i = 2:length(image_grid)
    restrictFlag = image_grid(i).restrict;
    imFixed = array_restrict(imFixed,restrictFlag);
    nanrng = nanbox(imFixed);
    if any(diff(nanrng,1,2) < options.min_pixels)
      i = i-1; % for setting gap_data below
      break
    end
    image_grid(i).imFixed = imFixed;
    mask = rmo_restrictmask(mask,restrictFlag);
    image_grid(i).mask = mask;
    g0c = register_g0(image_grid(i).sz,imclass);
    image_grid(i).g0 = cat(n_dims+1,g0c{:});
  end
  
  options.image_grid = image_grid;
  if strcmp(options.gap_data_mode,'auto')
    options.gap_data = length(image_grid)-i;
  end
  
function maskr = rmo_restrictmask(mask,restrictFlag)
  maskd = double(mask);
  maskdr = array_restrict(maskd,restrictFlag);
  maskr = (maskdr > 1-10*eps);
  maskr = imerode(maskr,ones(repmat(3,1,ndims(mask))));
