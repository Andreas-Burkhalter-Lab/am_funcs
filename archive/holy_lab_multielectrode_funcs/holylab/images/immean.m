function imout = immean(varargin)
% IMMEAN takes the mean over a set of images
% IMSUM: compute the mean value of images
% Syntax:
%   imout = immean(im1,im2,...)
% where
%   im1, im2, ... are the input images
%   imout is the resultant average. This will be of type 'single' no
%     matter what format the inputs are in.
  
  if (nargin < 1)
    imout = [];
    return;
  end
  imout = single(varargin{1});
  for i = 2:length(varargin)
    imout = imout + single(varargin{i});
  end
  imout = imout/length(varargin);
  
  
  