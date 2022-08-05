function varargout = imqinterp(g,varargin)
% IMQINTERP: quadratic interpolation and gradients for images
%
% This function is designed to facilitate optimization problems that
% involve of functions of images that have to be evaluated off-grid
% (i.e., with sub-pixel accuracy).  It uses a piecewise quadratic
% interpolation scheme (uniform quadratic B-spline resampling),
% resulting in a function that is continuous in both value and gradient.
%
% In addition to returning interpolated images, this function can also
% compute the gradient of the interpolated function with respect to
% changes in the evaluation location.  This gradient is frequently useful
% in optimization problems.
% 
% A common use-case for this function is image registration.
%
% Syntax:
%   [im1o,im2o,...] = imqinterp(g,im1,im2,...)
% Evaluates im1, im2, ... at the positions listed in g. Here are the
% details (a simple example is shown below):
%   g is an array of size [sz_out n_dims_in], where sz_out is a vector
%     specifying the size of the output images (subject to caveats
%     described below), and n_dims_in is the # of spatial dimensions in
%     the input images.
%   The input images, im1 and higher, have at least n_dims_in dimensions,
%     but may have additional dimensions.  For example, if im1 is an RGB
%     image it would be of size [sz1 3], where the last coordinate
%     corresponds to the color channel and length(sz1) == n_dims_in.
%     The number of dimensions might vary among the different
%     input images.  However, the size of the spatial part must be
%     uniform across input images.
%   The output images, im1o and higher, are arrays of (at least) size
%     sz_out, but will also have additional dimensions corresponding to any
%     additional dimensions of the inputs.
%
% NOTE: if your input images are one-dimensional, g must be arranged as a
% column vector.
%
% Examples:
%   1. im1 and im2 are of size [nx ny nz] and g is of size [nx ny nz 3].
%      This corresponds to a case in which two grayscale 3-dimensional
%      images are being evaluated off-grid at "full size."
%      If im2 were instead an RGB image, it would have size
%      [nx ny nz 3], and each color channel would be interpolated
%      independently for the RGB output image im2o.
%   2. Given a two-dimensional input image im1 of size [nx ny], one could
%      interpolate just a subregion if g is of size [nxo nyo 2], where
%      nxo <= nx and nyo <= ny.
%   3. If im1 is two dimensional and g is of size [n 2], the interpolated
%      image would simply be evaluated at a list of locations; in this
%      case, the output image would be a column vector (or of size [n 3]
%      if im1 were an RGB image).
%   4. One could compute a two-dimensional "slice" of a three-dimensional
%      image by having im1 be of size [nx ny nz] and g be of size [n1 n2 3].  
%      The output two-dimensional image would be of size [n1 n2].
%      The slice could be at any angle, not necessarily aligned with a
%      coordinate axis, and would not even need to be planar.
% 
% One can also get the gradients with respect to the warping function g,
% by requesting additional outputs:
%   [im1o,im2o,...,im1grad,im2grad,...] = imqinterp(g,im1,im2,...)
% where additionally im1grad, im2grad, ... are arrays of at least size
%   [sz_out n_dims_in], but may have additional dimensions if the input
%   images do.  The additional dimensions come _after_ the gradient
%   direction dimension.
%
% Example:
%   im = double(imread('cameraman.tif'));
%   % Make g, the deformation. This will be a random deformation,
%   % so the results will be very weird!
%   [X1,X2] = ndgrid(1:size(im,1),1:size(im,2));
%   g = X1; g(:,:,2) = X2;
%   g = g + randn(size(g));
%   % Do the interpolation
%   imi = imqinterp(g,im);
%   figure; subplot(1,2,1); imshowsc(im); subplot(1,2,2); imshowsc(imi)
%
% See also: IMINTERP, IMTRANSFORM.
  
% Copyright 2010 by Timothy E. Holy and Zhongsheng Guo
  
% Implemented as a MEX file.
