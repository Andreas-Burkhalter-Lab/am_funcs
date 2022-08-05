function options = register_phasecorr_initialize(im,options)
% register_phasecorr_initialize: set parameters for phasecorr registration
% Syntax:
%   options = register_phasecorr_initialize(im)
%   options = register_phasecorr_initialize(im,options)
% where
%   im is an image. im may be two-dimensional or higher, and it can
%     optionally have multiple values for each pixel (e.g., color
%     channels).
%   options is a structure which may have the following fields:
%     mode: 'grayscale', 'rgb' (or 'single-valued'/'multi-valued') used
%       to indicate whether the image has a single value per pixel or
%       multiple values per pixel. If multi-valued, the last dimension of
%       the array contains the values, so that (for example) a
%       two-dimensional rgb image would have size m-by-n-by-3.  Default:
%       grayscale.
%     overlap_fraction (default 0.33): the phase correlation algorithm
%       works by cutting out "chunks" of the image, one per grid element of
%       the deformation u. It then finds the best translation (using the
%       normalized phase correlation) to align each chunk to the
%       corresponding portion of the fixed (non-moving) image. If
%       overlap=0, the chunks are entirely independent; if overlap>0, then
%       each chunk will be larger, including regions of overlap (i.e.,
%       redundancy) with neighboring chunks.
%       0.33 means that the overlap is 33% of the gap between the centers
%       of adjacent chunks.
%     smooth: if present, smooths the normalized correlation before
%       searching for the peak of the correlation. This can be supplied in
%       one of two ways: first, you can use a function handle of the form
%           rs = func(r)
%       where r is the normalized correlation (i.e., an array) and rs is
%       the smoothed correlation. One example would be
%           func = @(r) imfilter(r,h)
%       where h is a smoothing kernel.
%       Alternatively, you can specify this as a vector "sigma" for
%       gaussian smoothing. For example, for a two-dimensional image you
%       might specify this as [1 1], and for a three-dimensional stack you
%       might consider [1 1 0] (no smoothing in the third coordinate). For
%       sigmas that are small, you may get better performance from
%       imfilter.
%       By default, this field is not set.
%     lambda (default 0): set to positive if you want to penalize changes
%       in volume under deformation.  This is in units of (pixel value)^2,
%       so the mean square pixel value of your image might give you a
%       reasonable indication of what values of lambda correspond to a
%       large or small penalty.
%     pixel_spacing and/or u_schedule OR pyramid: use these if you want to use
%       multiresolution/multigrid techniques, e.g., to start with a coarse
%       deformation (e.g., a global shift) and then successively refine it
%       to more and more local detail. These fields control the schedule of
%       refinement of the deformation to finer grids.
%          If you specify u_schedule, you explicitly give the size of the grid
%       at each stage of refrinement.  For example (assuming two dimensional
%       images),
%              options.u_schedule = [1 1; 2 2; 3 3; 5 5; 9 9; 17 17];
%       would cause it to start with a rigid shift, then on the next call
%       progress to a 2-by-2 grid, then a 3-by-3 grid, a 5-by-5 grid, and
%       so on.
%          If you specify pixel_spacing (as a vector giving the displacement
%       in physical units between pixels along each coordinate), then the
%       refinement of the u-grid will progress in a way that eventually
%       causes the chunks cut from the image will be designed to be of
%       approximately equal physical size along each dimension.
%          Finally, a field "pyramid" specifies the schedule as a structure
%       array of the form of the output array_restrict_schedule. This form
%       has a significant performance advantage (it allows use of
%       array_prolong when creating the warped image), and should therefore
%       be used whenever possible.
%          u_schedule takes precedence over pixel_spacing, until it is
%       exhausted, after which pixel_spacing will control additional
%       refinements (if you keep asking for them).
%          If you specify none of these, the default is to increase to 2*sz-1,
%       where sz is the previous grid size.
%          These parameters are used by register_phasecorr_prolong.
%
% See also: register_phasecorr_prolong, array_restrict_schedule.

% Copyright 2010 by Timothy E. Holy

  %% Initialize
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'mode','grayscale','overlap_fraction',0.33);
  
  % Determine the dimensionality
  switch lower(options.mode(1:3))
    case {'gra','sin'}
      n_dims = ndims(im);
      sz = size(im);
      sz_spatial = sz;
      n_values = 1;
    case {'mul','rgb'}
      sz = size(im);
      sz_spatial = sz(1:end-1);
      n_values = sz(end);
      n_dims = length(sz_spatial);
  end
  % Correct for one-dimensional "images"
  if (sz_spatial(end) == 1)
    sz_spatial = sz_spatial(1:end-1);
    n_dims = n_dims-1;
  end
  
  % If necessary, provide the identity map
  if ~isfield(options,'g0')
    cl = class(im);
    if isempty(strmatch(cl,{'single','double'},'exact'))
      cl = 'single';
    end
    options.g0 = register_g0(sz_spatial,cl);
  end
  
  options.n_dims = n_dims;
  options.sz_spatial = sz_spatial;
  options.n_values = n_values;
  options.class = class(im);
  
  if ~isfield(options,'pyramid')
    % If we don't have a pyramid, set up an initial u-schedule to go from
    % global rigid translations to our first "real" grid, 2-by-2-by-...
    % From that point onward, either pixel_spacing or the default 2*n-1
    % behavior will govern refinement.
    options = default(options,'u_schedule',repmat([1;2],1,options.n_dims));
  end
    