function imshowrg(im1,im2,clim)
% IMSHOWRG: display an overlay of two images in red/green channels
% This is useful for testing image registration.
%
% Syntax:
%   imshowrg(im1,im2)
%   imshowrg(im1,im2,clim)
% where
%   im1 and im2 are grayscale images of the same size;
%   clim (optional) are the color limits, a 2-vector giving the values
%     that map to 0 and 1, respectively (default: taken from max & min
%     values).
%
% See also: IMSHOW, IMSHOWSC.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 2)
    im2 = [];
  end
  if (ndims(im1) > 2 || ndims(im2) > 2)
    error('Works only with grayscale images');
  end
  if (nargin < 3)
    if ~isempty(im2)
      clim = [min(double(min(im1(:))),double(min(im2(:)))) max(double(max(im1(:))),double(max(im2(:))))];
    else
      clim = double([min(im1(:)) max(im1(:))]);
    end
  end
  if (~isequal(size(im1),size(im2)) && ~isempty(im2))
    error('Works only when the images are of the same size');
  end
  % May need to magnify, by copying pixels
  sz = size(im1);
  screen_size = get(0,'ScreenSize');
  target_size = min(screen_size(3:4));
  if (max(sz) < target_size/4)
    n_rep = round(target_size/max(sz));
    for dimIndex = 1:2
      tmp = repmat(1:sz(dimIndex),n_rep,1);
      indx{dimIndex} = tmp(:);
    end
    im1 = im1(indx{:});
    if ~isempty(im2)
      im2 = im2(indx{:});
    end
  end
  if ~isempty(im2)
    im = zeros([size(im1),3],'single');
    im(:,:,1) = imshowrg_scale(im1,clim);
    im(:,:,2) = imshowrg_scale(im2,clim);
  else
    im = imshowrg_scale(im1,clim);
  end
  imshow(im);
  
  
function imsc = imshowrg_scale(im,clim)
  imsc = (single(im) - clim(1))/diff(clim);
  imsc(imsc(:) < 0) = 0;
  imsc(imsc(:) > 1) = 1;
  