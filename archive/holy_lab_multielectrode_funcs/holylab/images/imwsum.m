function imo = imwsum(weight,varargin)
% IMWSUM: calculate the weighted sum of images
%
% If w(i) is a scalar and im(i) is an image, this function computes
%    sum_i(w(i)*im(i))
%
% Syntax:
%   imout = imwsum(w,A,B,C,...)
% where
%   weight is the vector of weights;
%   A,B,C,... are the images (note length(weight) == # images)
% and
%   imout is the result.
%
% See also: IMMEAN.
  
% Copyright 2005 by Timothy E. Holy
  
  nimages = length(varargin);
  if (nimages ~= length(weight))
    error(['The number of weights does not agree with the number of ' ...
           'images']);
  end
  if (nimages == 0)
    imo = [];
    return;
  end
  
  imo = weight(1)*single(varargin{1});
  for i = 2:nimages
    imo = imo + weight(i)*single(varargin{i});
  end
  