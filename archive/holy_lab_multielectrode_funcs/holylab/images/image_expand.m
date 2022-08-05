function [im_out expand_matrices] = image_expand(im_in, szout)
% image_expand takes a small image and expands the exact image to a larger
% size (for any reason you might like)
% Syntax im_out = image_expand(im_in, szout)
%
% Where:
%   im_in is the image you wish to expand
%   szout is the dimensionality to which you want to expand the image
% and
%   im_out is the expanded image
%   expand_matrices is a cell array containing the width and height
%                   expansion matrices, respectively
%   e.g., im_out = ((im_in*expand_matrices{1})'*expand_matrices{2})'
%
% Copyright 2010 Julian P. Meeks

%% Identify the dimensions of the input image
sz = size(im_in);
ncolin = sz(2);
nrowin = sz(1);

%% construct the width expansion matrix
% size(expand_matrices{1},2) must be equal to ncolin
ncolout = szout(2);
expand_matrices{1} = zeros([ncolin ncolout]);
% as closely as possible, create an even expansion matrix
step = ncolout/ncolin;
for i = 1:ncolin
    expand_matrices{1}(i,1+round(step*(i-1)):round(step*(i))) = 1; 
end
temp = im_in*expand_matrices{1}; temp = temp';
%% construct the height expansion matrix
nrowout = szout(1);
expand_matrices{2} = zeros([nrowin nrowout]);
step = nrowout/nrowin;
for i = 1:nrowin
    expand_matrices{2}(i,1+round(step*(i-1)):round(step*(i))) = 1; 
end
temp = temp*expand_matrices{2};
im_out = temp';
end