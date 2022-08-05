function [imf,hfft] = imfilter_fft(im,h0,center)
% IMFILTER_FFT: filter a multidimensional image using FFTs
% Syntax:
%   imf = imfilter_fft(im,h)
%   imf = imfilter_fft(im,h,center)
%   [imf,hfft] = imfilter_fft(...)
% where
%   im is the input image/stack/multidimensonal object;
%   h is the filter you want to use;
%   center is the coordinate (in indices) of the center position within
%     the filter (default: the center of h);
% and
%   imf is the filtered image;
%   hfft is the fourier transform of the filter (can be used for
%     deconvolution).
%
% Data type support: singles and doubles.
%
% This is to be preferred over IMFILTER when the size of the filter is
% large.
%
% See also: IMFILTER.
  
% Copyright 2006 by Timothy E. Holy
  
  nd = ndims(h0);
  if (nargin < 3)
    szh = size(h0);
    center = ceil(szh/2);
  end
  h = zeros(size(im),'single');
  x = cell(1,nd);
  % Create a wrap-around filter, so that it's centered on the center
  for i = 1:nd
    x{i} = [[1:center(i)-1]+size(im,i)-center(i)+1 (center(i):size(h0,i))-center(i)+1];
  end
  h(x{:}) = h0;
  % Take fourier transforms
  imfft = fftn(im);
  hfft = fftn(h);
  imf = ifftn(imfft.*hfft);
  