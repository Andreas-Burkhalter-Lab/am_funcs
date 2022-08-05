function err = im_gaussians_err(im,model,im_coords)
% IM_GAUSSIANS_ERR: compute the fitting error between objects and gaussians
% Syntax:
%   err = im_gaussians_err(im,model,im_coords)
% where
%   im is the actual image (i.e., the data);
%   model and im_coords define "patches" of the modelled image, as output
%     by GAUSSIANS2IM.
% and
%   err is a 1-by-n_objects vector, containing the fitting error for each
%     object. Note that, unlike the error returned by IM2GAUSSIANS, this
%     error measure includes all sufficiently-nonzero model pixels, and
%     thus may spill over the boundaries of the data object.
%
% See also:  GAUSSIANS2IM, IM2GAUSSIANS.

% Copyright 2006 by Timothy E. Holy

  n_objects = length(model);
  
  err = zeros(1,n_objects);
  for objIndex = 1:n_objects
    err_im = im(im_coords{objIndex}) - model{objIndex};
    err(objIndex) = sum(err_im(:).^2);
  end