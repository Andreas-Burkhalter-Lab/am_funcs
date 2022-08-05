function [dx,im2shift,err] = register_overlap_multires(im1,im2,options)
% REGISTER_TRANSLATION_MULTIRES: register by translation only (rigid)
% This uses a multi-resolution approach to efficiently and robustly
% calculate the shift to bring two images into optimal alignment.
%
% Syntax:
%   [dx,im2s] = register_translation_multires(im1,im2);
%   [dx,im2s,err] = register_translation_multires(im1,im2,options)
% where
%   im1 is the reference image (best if of type single or double)
%   im2 is the image to be shifted
%   options is a structure which may have the following fields:
%     dx_start (default zeros): starting guess for the appropriate shift
%       (see below for more detail). You probably only need to supply this
%       if the shift is fairly large on the scale of coarse image features,
%       and the algorithm has a hard time finding an appropriate solution.
%     pixel_spacing (default ones): the grid spacing, in physical units,
%       along each axis of the image. If your images are sampled with
%       different resolution along different axes, specifying the spacing
%       will allow the algorithm to choose a multiresolution pyramid that
%       tends to equalize the resolution as one proceed to coarser scales.
%       It's highly recommended you make use of this parameter if your
%       images do have different sampling along different axes.
%     min_pixels (default 15): as part of the multiresolution strategy,
%       images will be coarsened, but will not be made smaller than
%       min_pixels along each coordinate
%     n_vcycles (default 2): the number of V-cycles in descending &
%       ascending the grid hierarchy
%     display (default false): if true, causes fitting progress to be shown
%       as text
% and
%   dx is the amount of shift applied to the coordinates of im2 to bring it
%     in optimal alignment with im1
%   im2s is the shifted version of im2
%   err is a vector, err(i) is the mean square mismatch going into the ith
%     vcycle. err(end) contains the error at the end of the last v-cycle.

% Copyright 2007 by Timothy E. Holy, adapted for montaging by Julian P.
% Meeks

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'min_pixels',15);
  options.covariant = false;
  options = default(options,'n_vcycles',2);
  
  [overlap_img1, overlap_img2] = img_overlap(im1, im2, options.dx_start);
  options = register_multigrid_options(overlap_img1,options);
  options = default(options,'dx_start',zeros(1,options.n_dims));
  
  options.dx_start = [0 0];
  dx = options.dx_start;
  
  err = imcompare(overlap_img1,{overlap_img2});
  for iter = 1:options.n_vcycles
    [dx,im2shift,err(iter+1)] = register_translation_vcycle(dx,overlap_img2, ...
						  options);
    options.im2shift = im2shift;
  end
  
  