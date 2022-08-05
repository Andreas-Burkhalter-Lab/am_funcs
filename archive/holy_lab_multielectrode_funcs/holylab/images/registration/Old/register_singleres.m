function [g,img] = register_singleres(im1,im2,lambda,options)
% REGISTER_singleres: multi-resolution non-rigid image registration
% This is the function to use for a "new" registration problem, where you
% don't start with a good initial guess for the deformation.  See
% REGISTER_GRADDESCENT when you want to tweak an existing solution.
%
% Syntax:
%   [g,img] = register_singleres(im1,im2)
%   [g,img] = register_singleres(im1,im2,lambda,options)
% where
%   im1 is the "fixed" image;
%   im2 is the "moving" image;
%   lambda (default 0) is the regularization coefficient (you probably want
%     to try to get things to work with lambda = 0; this parameter usually
%     doesn't help);
%   options is a structure which may have the following fields:
%     g_decimate (default []) specifies a different schedule of image
%       reduction.  You only need specify the schedule for as many
%       reductions as you care to use; unspecified reductions will be by
%       threefold until the image is reduced to a single point.  This
%       parameter is specified as an nsteps-by-ndims matrix, where each row
%       specifies a decimation factor in each coordinate---see IMREDUCE.
% and
%   img is the registered im2;
%   g is the deformation (of the format described in REGISTER_G0; it can
%     be plotted with REGISTER_PLOTG; and it can be expanded with
%     REGISTER_EXPANDG).
%
% A complete example:
%   im1 = single(imread('spine.tif')); figure; imshowrg(im1)
%   im2 = imrotate(im1,5,'bilinear','crop'); figure; imshowrg(im1,im2)
%   [g,img] = register_singleres(im1,im2);
%   figure; imshowrg(im1,img)
%   register_plotg(g)
%
% See also: REGISTER_GRADDESCENT, REGISTER_G0, REGISTER_PLOTG,
%   REGISTER_EXPANDG, IMREDUCE.
  
% Copyright 2006 by Timothy E. Holy
  
  n_dims = ndims(im1);
  if (nargin < 3)
    lambda = 0;
  end
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'g_decimate')
    options.g_decimate = zeros(0,n_dims);   % An empty schedule
  end
  
  n_decimation_steps = size(options.g_decimate,1);
  psi1 = sqrt(single(im1));  % Everything must be of type single!
  psi2 = sqrt(single(im2));
  psi1tmp = psi1;
  for decIndex = 1:n_decimation_steps
    psi1tmp = imreduce(psi1tmp,options.g_decimate(decIndex,:)); % slow, but certain to track modifications in imreduce
  end
  if (ndims(psi1tmp) ~= ndims(psi1))
    error('register::wrongdim',...
      ['The decimation schedule changed the dimensionality of the image!\n' ...
      'Alter the decimation schedule so that the final result has the same dimensionality.']);
  end
  g_sz = size(psi1tmp);
  
  % Set up default stepsize
  mu = 0.1;
  % Set up default deformation (identity map, i.e., no deformation).
  % This guess might need to be modified if the two images are of
  % different sizes (g0 should then "point" to the "cropped"
  % region). Allow g0 to be supplied as an option?
  g = register_g0(g_sz);
  % Do the registration
  options.sqrt = true;
  options.output_size = size(psi1);  % size to expand g to
  % Optimize the deformation
  g = register_nonrigid(psi1,psi2,g,lambda,mu,options);
  if (nargout > 1)
    % Calculate the output image
    g_hires = register_expandg(g,size(im1)); % scale up deformation
    img = register_warp(im2,g_hires); % Calculate warped image
  end
