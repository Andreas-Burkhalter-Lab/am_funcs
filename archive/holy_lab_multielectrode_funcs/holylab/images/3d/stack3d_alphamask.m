function alpha = stack3d_alphamask(im,alpha_fac,alpha_clim)
% STACK3D_ALPHAMASK: pixelwise adjustment of alpha scaling in 3d volumetric imaging 
%
% The main use case for this function is to highlight a particular
% structure or region in an image stack, by creating an alpha_func for
% stack3d. You would use it like this:
%    alpha_func = @(im,dfof) stack3d_alphamask(im,alpha_fac)
% or
%    alpha_func = @(im,dfof) stack3d_alphamask(im,alpha_fac,alpha_clim)
% (The default value of alpha_clim is [0 1])
%
% Here, alpha_fac is an array of the size of im, with values between 0 and
% 1, that will be used to create the final transparency value like this:
%     alpha = alpha_fac .* alpha0
% where alpha0 is computed from the normalized intensity in the "normal" way,
%     alpha0 = (im - alpha_clim(1))/diff(alpha_clim)
%
% The alpha_func would be supplied as an option to stack3d.
%
% See also: STACK3D.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 3)
    alpha_clim = [0 1];
  end
  alpha = (im - alpha_clim(1))/diff(alpha_clim);
  alpha(alpha<0) = 0;
  alpha(alpha>1) = 1;
  alpha = alpha_fac.*alpha;
end
