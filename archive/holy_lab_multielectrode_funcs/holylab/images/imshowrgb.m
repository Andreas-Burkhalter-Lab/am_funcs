function imo = imshowrgb(varargin)
% IMSHOWRG: display an overlay of two images in red/green channels
% This is useful for testing image registration.
%
% Syntax (examples):
%   imshowrgb(im1,im2)
%   imshowrgb(im1,im2,im3,clim)
%   im = imshowrgb(...)
% where
%   im1 and im2 (etc.) are grayscale images of the same size;
%   clim (optional) are the color limits, a 2-vector giving the values
%     that map to 0 and 1, respectively (default: taken from max & min
%     values).
% and
%   the output 'im' is the rgb image.
%
% See also: IMSHOW, IMSHOWSC.
  
% Copyright 2006 & 2009 by Timothy E. Holy
  
  n_images = length(varargin);
  if ~isequal(size(varargin{1}),size(varargin{end}))
    n_images = n_images-1;
    clim = varargin{end};
  else
    % We have to compute clim, set to be the min & max
    clim = [Inf -Inf];
    for i = 1:n_images
      clim(1) = min(clim(1),min(varargin{i}(:)));
      clim(2) = max(clim(2),max(varargin{i}(:)));
    end
  end
  if (n_images > 3)
    error('Cannot supply more than 3 images');
  end
  im = zeros([size(varargin{1}),3],'single');
  for i = 1:n_images
    im(:,:,i) = imshowrg_scale(varargin{i},clim);
  end
  image(im);
  set(gca,'DataAspectRatio',[1 1 1],'TickDir','out');
  if (nargout > 0)
    imo = im;
  end
  
  
function imsc = imshowrg_scale(im,clim)
  imsc = (single(im) - clim(1))/diff(clim);
  imsc(imsc(:) < 0) = 0;
  imsc(imsc(:) > 1) = 1;
  