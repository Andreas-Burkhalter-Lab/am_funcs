function masko = imdilatend(mask)
% IMDIALEND: multidimensional single-pixel dilation of images
% Matlab's imdilate dilates only in the x-y plane, and replicates this in
% the z-dimension. This function performs dilation in all coordinates.
% However, it dilates by only a single pixel.
%
% Syntax:
%   masko = imdilatend(mask)
% where
%   mask is the input (logical) array;
% and
%   masko is the output logical array, where a given pixel will be 1 if:
%     (1) that pixel in mask is 1; or
%     (2) any nearest-neighbor of that pixel is 1 in mask.
% Nearest-neighbors differ from the given pixel in only one coordinate,
% with a difference of 1. This therefore corresponds to the "diamond"
% pattern in strel/imdilate.
%
% See also: IMDILATE, STREL.

% Copyright 2006 by Timothy E. Holy

n_dims = ndims(mask);
sz = size(mask);
colons = {':'};
colons = repmat(colons,[1 n_dims]);
masko = mask;
for dimIndex = 1:n_dims
    cIndex1 = colons;
    cIndex1{dimIndex} = 1:sz(dimIndex)-1;
    cIndex2 = colons;
    cIndex2{dimIndex} = 2:sz(dimIndex);
    masko(cIndex1{:}) = masko(cIndex1{:}) | mask(cIndex2{:});
    masko(cIndex2{:}) = masko(cIndex2{:}) | mask(cIndex1{:});
end
