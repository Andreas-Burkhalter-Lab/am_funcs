function cm = imrange_to_colormap(range,cmsize)
% IMRANGE_TO_COLORMAP: set up a linear colormap given a range of intensities
% Syntax:
%   cm = imrange_to_colormap(range,cmsize)
% where
%   range is a 2-vector, [pixmin pixmax], giving the range of pixel
%     intensities over which the colormap is to be linear (and saturated
%     beyond)
%   cmsize is a scalar, giving the size of the colormap
% and
%   cm is a colormap, a cmsize-by-3 matrix in the range from 0 to 1
cm = zeros(cmsize,3);
range = double(range);
for i = range(1)+1:range(2)+1
  cm(i,:) = (i-range(1))/diff(range);
end
cm(range(2)+1:end,:) = 1;
