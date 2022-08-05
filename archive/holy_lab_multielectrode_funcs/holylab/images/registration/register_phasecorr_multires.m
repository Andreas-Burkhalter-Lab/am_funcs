function [u,imw,err] = register_phasecorr_multires(fixed,moving,options,multires_params)
% register_phasecorr_multires: iterative refinement of deformation until improvement stops
% Syntax:
%   [u,imw,err] = register_phasecorr_multires(fixed,moving,options,multires_params)
% where
%   fixed, moving, and options are as described in
%     register_phasecorr_improve
%   multires_params is a structure with the following fields:
%     desired_ulevel (required): the index into options.pyramid that
%       specifies the desired resolution of u.
%     coarsest_ulevel (default length(options.pyramid)): the coarsest u-resolution
%       at which this function is allowed to terminate. Early termination
%       occurs if refining the u-grid does not yield improvement, unless
%       this coarsest resolution has not yet been reached. To prevent early
%       termination, set coarsest_ulevel = desired_ulevel.
%     display (default false): if true, displays information on the
%       progress.
% and
%   u is the deformation. It is possible that it will be of lower resultion
%     than indicated by desired_ulevel, if the error stopped decreasing.
%   imw is the warped moving image
%   err is a vector of pixel mismatches at each scale of the pyramid.
%
% See also: register_phasecorr_refine.

% Copyright 2011 by Timothy E. Holy

  multires_params = default(multires_params,'coarsest_ulevel',length(options.pyramid),...
    'display',false);
  if isfield(options,'shift')
    imw = image_shift(moving,options.shift);
  else
    imw = moving;
  end
  err = nanmean((fixed(:)-imw(:)).^2);
  u = [];
  ulevel = length(options.pyramid)+1;
  while (ulevel >= multires_params.desired_ulevel)
    % Save the previous solution, in case the next is not an improvement
    uold = u;
    imwold = imw;
    % Refine the u-grid
    u = register_phasecorr_refine(u,fixed,moving,options);
    % Determine whether this is an improvement
    imw = register_phasecorr_warp(u,moving,options);
    thiserr = nanmean((fixed(:)-imw(:)).^2);
    if (thiserr < err(end) || ulevel >= multires_params.coarsest_ulevel)
      err(end+1) = thiserr; %#ok<AGROW>
    else
      if multires_params.display
        fprintf('Early termination, err = %g\n',thiserr);
      end
      u = uold;
      imw = imwold;
      break
    end
    if multires_params.display
      fprintf('u grid:');
      usz = size(u{1});
      usz(end+1:options.n_dims) = 1;
      fprintf(' %d',usz);
      fprintf(', error: %g\n',err(end));
    end
    ulevel = ulevel-1;
  end
  
  