function f = deconv_gaussian(width,noise)
% DECONV_GAUSSIAN: create a deconvolution filter for gaussian blur
% 
% Syntax:
%   f = deconv_gaussian(width,noise)
% where
%   width is the standard deviation of the gaussian, in # of pixels;
%   noise is the combination of readout + shot noise, in digital numbers
%     (i.e., pixel values) (default 0)
% and
%   f is the optimal deconvolution filter
  
% Copyright 2006 by Timothy E. Holy
  
  half_size = ceil(3*width);
  if (nargin < 2)
    noise = 0;
  end
  
  x = linspace(-half_size,half_size,2*half_size+1);
  m = exp(-x.^2/(2*width^2));  % The "measured" intensity profile
  m = m/sum(m);  % normalize
  
  f = deconv_filter_create(m,noise);
  