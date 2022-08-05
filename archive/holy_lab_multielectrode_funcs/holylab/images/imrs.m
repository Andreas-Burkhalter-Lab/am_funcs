function imout = imrs(im,n_pixels,flag)
% IMRS: image resample at lower resolution
% This is useful for "faithfully" representing large images (or matrices)
% under conditions where the number of pixels in the display is
% substantially smaller than the image or matrix.
%
% Syntax:
%   imout = imrs(im,n_pixels)
% where
%   im is the raw image or matrix;
%   n_pixels is a 2-vector containing the desired number of pixels along
%     each coordinate.  The image/matrix will be downsized (with
%     appropriate anti-aliasing) to be no bigger than a factor of 2 along
%     each of these dimensions.
%
%   imout = imrs(im,n_pixels,'lock')
%   locks the aspect ratio, so that the image/matrix will only be
%   downsized along one coordinate if it is also to be downsized in the
%   other.
%
% See also: ARRAY_RESTRICT.
  
% Copyright 2009 by Timothy E. Holy
  
  if (nargin < 2)
    n_pixels = [256 256];
  end
  lock = false;
  if (nargin > 2)
    if ischar(flag) && strcmp(flag,'lock')
      lock = true;
    end
  end
  
  imout = im;
  resflag = size(imout) > 2*n_pixels;
  downsize = any(resflag);
  if lock
    downsize = all(resflag);
  end
  while downsize
    imout = array_restrict(imout,resflag);
    resflag = size(imout) > 2*n_pixels;
    downsize = any(resflag);
    if lock
      downsize = all(resflag);
    end
  end
end
