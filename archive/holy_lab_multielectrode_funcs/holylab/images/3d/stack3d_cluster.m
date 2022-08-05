function h = stack3d_cluster(im,clust,options)
% stack3d_cluster: highlight a set of pixels in a cluster
%
% Syntax:
%   h = stack3d_cluster(im,clust,options)
% where
%   im is a 3-dimensional stack
%   clust is a logical array of the same size as im, true for the pixels in
%     the cluster
%   options is a structure which may have the following fields:
%     nonclust_fac (default 0.1): the factor by which pixels outside the
%       cluster are "dimmed" (made more transparent). Smaller values make
%       them more transparent
%     Plus, any fields that stack3d (and therefore vol3d) accept (e.g.,
%     pixel_spacing, parent, etc.)
%
% There is a fully-worked out visualization example in the help for
% spin_flip_timer.
%
% See also: STACK3D, VOL3D, SPIN_FLIP_TIMER.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  if ~isequal(size(im),size(clust))
    error('Image and cluster are not of the same size');
  end
  options = default(options,'alpha_clim',[0 1],'nonclust_fac',0.1);
  afac1 = ones(size(im));
  afac2 = options.nonclust_fac*ones(size(im)); afac2(clust) = 1;
  options.alpha_func = @(im,dfof) stack3d_alphamask(im,afac1,options.alpha_clim);
  h = stack3d(im,options);
  set(h.handles,'HandleVisibility','off','Visible','off');
  options.alpha_func = @(im,dfof) stack3d_alphamask(im,afac2,options.alpha_clim);
  h(2) = stack3d(im,options);
end
