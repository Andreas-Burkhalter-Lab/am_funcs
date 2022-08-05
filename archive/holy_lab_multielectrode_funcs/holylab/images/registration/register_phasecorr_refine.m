function [uout,imwold] = register_phasecorr_refine(u,fixed,moving,options)
% register_phasecorr_refine: change a coarse deformation to a finer one
%
% In performing registration, one of the best strategies is to progress
% from rigid deformations to warping deformations, with the warp grid
% iteratively becoming finer and finer.  This helps ensure that the
% solutions are as "smooth" as possible. This function is designed to
% make that process straightforward.
%
% Syntax:
%   [unew,imwold] = register_phasecorr_refine(uold,fixed,moving,options)
% where
%   uold is either [] (to start from scratch with a rigid shift) or a
%     deformation, for example from a previous call
%   fixed, moving, and options are as described in
%     register_phasecorr_improve
% and
%   unew is a cell array containing the new deformation at finer grid
%     sampling. Not only is it interpolated up from the old deformation,
%     but it is optimized by a single call to
%     register_phasecorr_improve.
%   imwold is as described in register_phasecorr_improve.
%
% Example:
%   % Step 1: create a warped image using a random deformation
%   im = double(imread('cameraman.tif'));
%   scale = 10;  % max # of pixels in shift
%   sz = size(im);
%   udef = randn([sz 2]);
%   udef = imfilter_gaussian(udef,[30 30 0]);  % smooth the random deformation
%   udef = udef * (scale / max(abs(udef(:)))); % set the scale
%   imw_full = register_multigrid_warp(im,udef);
%   rng = 15:240;
%   moving = imw_full(rng,rng);
%   fixed = im(rng,rng);
%   figure; imshowrgb(fixed,moving); title('Unregistered images')
%   % Step 2: perform the registration
%   options = register_phasecorr_initialize(fixed);
%   options.lambda = 1;  % penalize non-smooth deformations
%   options.smooth = [1 1];  % use pixel smoothing to calculate opt. shift
%   u = [];
%   u = register_phasecorr_refine(u,fixed,moving,options); % 1-by-1 (rigid)
%   u = register_phasecorr_refine(u,fixed,moving,options); % 2-by-2
%   u = register_phasecorr_refine(u,fixed,moving,options); % 3-by-3
%   u = register_phasecorr_refine(u,fixed,moving,options); % 5-by-5
%   u = register_phasecorr_refine(u,fixed,moving,options); % 9-by-9
%   u = register_phasecorr_refine(u,fixed,moving,options); % 17-by-17
%   imw = register_phasecorr_warp(u,moving,options);
%   figure; imshowrgb(fixed,imw); title(['Registered, grid size [' num2str(size(u{1})) ']'])
%
% See also: register_phasecorr_improve, register_phasecorr_initialize.

% Copyright 2010 by Timothy E. Holy
  
  [~,u_old] = register_phasecorr_prolong(u,options);
  [uout,imwold] = register_phasecorr_improve(u_old,fixed,moving,options);
  