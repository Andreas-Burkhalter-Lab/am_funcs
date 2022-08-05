function cm = graymap(limits)
% GRAYMAP: set up extended gray colormap
% Syntax:
%   cm = graymap(limits)
% where
%   limits is a 2-vector [blacklevel whitelevel]
%
% Using this colormap with 'image' is equivalent to imagesc(im,limits)
% and colormap(gray)
%
% See also: COLORMAP.
  
% Copyright 2005 by Timothy E. Holy
  
  cm = zeros(max(limits),3);
  cm(limits(1):limits(2),:) = repmat(linspace(0,1,diff(limits)+1)',1,3);
  