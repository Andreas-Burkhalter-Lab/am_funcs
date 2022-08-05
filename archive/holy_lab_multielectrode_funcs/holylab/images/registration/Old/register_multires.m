function [g,img,err] = register_multires(im1,im2,lambda,options)
% REGISTER_MULTIRES: multi-resolution non-rigid image registration
% This is the function to use for a "new" registration problem, where you
% don't start with a good initial guess for the deformation.  See
% REGISTER_NONRIGID when you want to tweak an existing solution.
%
% Syntax:
%   [g,img] = register_multires(im1,im2)
%   [g,img,err] = register_multires(im1,im2,lambda,options)
% where
%   im1 is the "fixed" image;
%   im2 is the "moving" image;
%   lambda (default 0) is the regularization coefficient (you probably want
%     to try to get things to work with lambda = 0; this parameter usually
%     doesn't help);
%   options is a structure which may have the following fields:
%     covariant (default true): perform covariant registration, in which
%       volume changes have a corresponding effect on intensity,
%       preserving total intensity over local regions.
%     pyramid_schedule (default []) specifies a different schedule of image
%       reduction.  You only need specify the schedule for as many
%       reductions as you care to use; unspecified reductions will be by
%       threefold until the image is reduced to a single point.  This
%       parameter is specified as an nsteps-by-ndims matrix, where each row
%       specifies a decimation factor in each coordinate---see IMREDUCE.
%     n_levels (default 3) is the # of pyramid steps by which the images
%       should be larger than the deformation grid.  Larger n_levels means
%       a coarser grid.  For example, if your image is 1024x1024, and you
%       are using pyramid reductions by a factor of 3, then
%       setting n_levels to be 3 would make the grid about 27x smaller, or
%       of size 38x38.
%     g0 allows you to specify a starting guess for g (see REGISTER_G0
%       for a description of the parametrization).
% and
%   img is the registered im2;
%   g is the deformation (of the format described in REGISTER_G0; it can
%     be plotted with REGISTER_PLOTG; and it can be expanded with
%     REGISTER_EXPANDG).
%   err is the discrepancy between the images.
%
% A complete example:
%   im1 = single(imread('spine.tif')); figure; imshowrg(im1)
%   im2 = imrotate(im1,5,'bilinear','crop'); figure; imshowrg(im1,im2)
%   [g,img] = register_multires(im1,im2);
%   figure; imshowrg(im1,img)
%   register_plotg(g)
%
% See also: REGISTER_NONRIGID, REGISTER_G0, REGISTER_PLOTG,
%   REGISTER_EXPANDG, IMREDUCE.
  
% Copyright 2006 by Timothy E. Holy
  
  n_dims = ndims(im1);
  if (nargin < 3)
    lambda = 0;
  end
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'n_levels')
    options.n_levels = 3;
  end
  if ~isfield(options,'pyramid_schedule')
    options.pyramid_schedule = zeros(0,n_dims);   % An empty schedule
  end
  if ~isfield(options,'covariant')
    options.covariant = true;
  end
  
  n_levels = options.n_levels;
  pyramid_schedule = options.pyramid_schedule;
  im1 = single(im1);  % Everything must be of type single!
  im2 = single(im2);
  indexPyr = 1;
  im1p{1} = im1;      % Set up pyramids; supplied images are the ...
  im2p{1} = im2;      % ... finest scale of analysis
  % Create the pyramid schedule to reduce the image all the way down to a
  % single point.
  % This might seem mildly wasteful, since we only have to go n_levels
  % above this point, but the amount of data is so small that it's not an
  % issue---and it really helps insure that the sizes of g will be correct.
  while (indexPyr <= size(pyramid_schedule,1) || ...
	 any(size(im1) > 1))
    if (indexPyr <= size(pyramid_schedule,1))
      % Use user-supplied step
      thisPyr = pyramid_schedule(indexPyr,:);
    else
      % Create a default step
      thisPyr = 3*ones(1,n_dims);
    end
    pyramid_schedule(indexPyr,:) = thisPyr;
    % Shrink the images
    im1 = imreduce(im1,thisPyr);
    im2 = imreduce(im2,thisPyr);
    im1p{end+1} = im1;
    im2p{end+1} = im2;
    indexPyr = indexPyr+1;
  end
  % Adjust the level separation for the user if needed
  n_levels_old = n_levels;
  while (ndims(im1p{end-n_levels}) ~= ndims(im1p{1}))
    n_levels = n_levels + 1;
  end
  if (n_levels ~= n_levels_old)
    warning('Number of levels between g & image data (n_levels) changed to preserve dimensionality');
  end
  % Set up default stepsize
  mu = 0.1;
  % Set up default deformation (identity map, i.e., no deformation).
  % This guess might need to be modified if the two images are of
  % different sizes (g0 should then "point" to the "cropped"
  % region). Allow g0 to be supplied as an option?
  have_g0_pyramid = false;
  if ~isfield(options,'g0')
    g = register_g0(ones(1,n_dims));
  else
    % The user supplied a guess for g. Create the pyramid, so we get g0's
    % appropriate for each level of analysis.
    g0_pyramid{1} = options.g0;
    for indexPyr = n_levels+1:size(pyramid_schedule,1)
      g0Index = length(g0_pyramid)+1;
      g0_pyramid{g0Index} = register_reduceg(g0_pyramid{g0Index-1},...
          pyramid_schedule(indexPyr,:));
    end
    g = g0_pyramid{end};
    have_g0_pyramid = true;
  end
	
  for indexPyr = length(im1p)-n_levels:-1:1
    if options.covariant
      psi1 = sqrt(im1p{indexPyr});  % use covariant registration
      psi2 = sqrt(im2p{indexPyr});  % "
      ops.sqrt = true;
    else
      psi1 = im1p{indexPyr};
      psi2 = im2p{indexPyr};
      ops.covariant = false;
      ops.sqrt = false;
    end
    ops.output_size = size(psi1);  % size to expand g to
    ops.g_decimate = pyramid_schedule(indexPyr:indexPyr-1+n_levels,:);
                        % How to shrink down the gradient
    % Optimize the deformation
    [g,psig,muOpt,err] = register_nonrigid(psi1,...
      psi2,...
      g,...
      lambda,...
      mu,...
      ops);
    % "Blow up" the deformation to the next scale of analysis
    if (indexPyr > 1)
      g = register_expandg(g,size(im1p{indexPyr+n_levels-1}));
      if have_g0_pyramid
          % Since we had a previous solution for g, use the additional
          % fine-scale structure between the two corresponding levels of g0
          % to correct our initial guess for the next scale of g
          g0_expand = register_expandg(g0_pyramid{indexPyr},...
              size(g0_pyramid{indexPyr-1}{1}));
          for dimIndex = 1:n_dims
              dg0 = g0_pyramid{indexPyr-1}{dimIndex} - g0_expand{dimIndex};
              g{dimIndex} = g{dimIndex} + 0.8.*dg0;
          end
      end
    end
  end
  
  if (nargout > 1)
    % Calculate the output image
    g_hires = register_expandg(g,size(im1p{1})); % scale up deformation
    img = register_warp(im2p{1},g_hires,ops); % Calculate warped image
  end
