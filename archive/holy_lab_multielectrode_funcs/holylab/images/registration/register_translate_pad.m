function impad = register_translate_pad(im)
% register_translate_pad: pads an image with NaNs in preparation for fourier transforms
%
% This function works with register_translate_nancorr to compute the
% optimum translation aligning two images.
%
% Syntax:
%   impad = register_translate_pad(im)
%
% See also:  register_translate_nancorr.

% Copyright 2011 by Timothy E. Holy

  sz = size(im);
  impad = nan(2*sz-1);
  coords = cell(1,length(sz));
  for i = 1:length(sz)
    coords{i} = 1:sz(i);
  end
  impad(coords{:}) = im;
  