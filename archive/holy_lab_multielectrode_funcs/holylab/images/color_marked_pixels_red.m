function imout = color_marked_pixels_red(im,mark)
% COLOR_MARKED_PIXELS_RED: display image in grayscale except chosen pixels are red
%
% Syntax:
%   imout = color_marked_pixels_red(im,mark)
% where
%   im is a 2-dimensional image
%   mark is a logical matrix of the same size as image, the entries that
%     are true will be colored red
% and
%   imout is a RGB array (the first two dimensions the same size as im)
%     that can be displayed with image or imagesc.

% Copyright 2009 by Timothy E. Holy, with input from Rachel O. Wong

  if (ndims(im) ~= 2)
    error('Input image must be two-dimensional')
  end
  imblanked = im;
  imblanked(mark) = 0;
  imout = repmat(imblanked,[1 1 3]);
  imout(:,:,1) = im;
end
