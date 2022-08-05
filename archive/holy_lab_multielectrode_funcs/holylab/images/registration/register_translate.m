function shift = register_translate(im1,im2,options)
% REGISTER_TRANSLATE: find the translation that best aligns images
%
%   shift = register_translate(im1,im2)
% where
%   im1 is the "fixed" image (may be of arbitrary dimensionality)
%   im2 is the "moving" image
% and
%   shift is a vector containing the "optimal" displacements.
%
%   shift = register_rigid(im1,im2,ops)
% lets you control behavior via some options. See array_findpeak for the
% available choices.
%
% Example:
%   im = imread('cameraman.tif');
%   x = 50:200;
%   dx = 15; dy = -20; f = 0.7;
%   im1 = im(x,x);
%   im2 = f*im(x+dx,x+dy) + (1-f) * im(x+dx+1,x+dy);
%   shift = register_translate(im1,im2)
%
% See also: REGISTER_RIGID, ARRAY_FINDPEAK.

% Copyright 2011 by Timothy E. Holy

  if isinteger(im1)
    im1 = single(im1);
  end
  if isinteger(im2)
    im2 = single(im2);
  end
  sz = size(im1);
  if ~isequal(sz,size(im2))
    error('Two images must be of equal size')
  end
  if (nargin < 3)
    options = struct;
  end
  
  % Compute the maximum of the phase correlation
  if ~isreal(im1)
    im1f = im1;  % the user already supplied this as the fft
  else
    im1f = fftn(im1);
  end
  im2f = fftn(im2);
  rf = im2f .* conj(im1f) ./ abs(im2f.*im1f);
  rf(isnan(rf)) = 0;
  r = ifftn(rf);
  shift = array_findpeak(fftshift(r),options);
