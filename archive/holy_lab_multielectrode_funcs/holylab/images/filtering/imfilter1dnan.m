% IMFILTER1DNAN: filter an image along one dimension, ignoring NaNs
%
% Syntax:
%   imout = imfilter1dnan(im,h,dim)
% where
%   im is the input image (may be multidimensional)
%   h is the 1-d filter
%   dim is the dimension to filter along
% and
%   imout is the filtered image.
%
% See also: IMFILTER, IMFILTER_GAUSSIAN.

% Copyright 2006 by Timothy E. Holy

% This is a MEX file.