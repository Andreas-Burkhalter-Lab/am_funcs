function masko = trim_mask_edges(mask,sz)
% trim_mask_edges: line a logical array with false "pixels"
% Syntax:
%   masko = trim_mask_edges(mask,sz)
% where
%   mask is a logical array
%   sz is a vector of length ndims(mask); sz(dimIndex) pixels are set to
%     false along the corresponding edge of masko.
% and
%   masko is the output mask with edges set to false.
%
% Example:
%   mask = imF > 2500;
%   mask = trim_mask_edges(mask,[15 15 3]);

% Copyright 2010 by Timothy E. Holy

  masko = mask;
  n_dims = ndims(mask);
  colons = repmat({':'},1,n_dims);
  for dimIndex = 1:n_dims
    cc = colons;
    cc{dimIndex} = 1:sz(dimIndex);
    masko(cc{:}) = false;
    cc{dimIndex} = cc{dimIndex} + size(mask,dimIndex)-sz(dimIndex);
    masko(cc{:}) = false;
  end