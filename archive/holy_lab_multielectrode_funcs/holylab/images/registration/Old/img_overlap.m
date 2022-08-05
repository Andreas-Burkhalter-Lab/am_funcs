function [overlap_img1, overlap_img2] = img_overlap(im1, im2, dx)
% img_overlap returns only the overlapping regions of two images
% intended for use with register_overlap_multires
%
% In progress:

% Copyright 2007 Julian P Meeks

% Set y-range of images
if dx(1) >= 0
    yrange1 = [1+dx(1) size(im1, 1)];
    yrange2 = [1 size(im2, 1)-dx(1)];
else
    yrange1 = [1-dx(1) size(im1,1)];
    yrange2 = [1 size(im2, 1)+dx(1)];
end

% Set x-range of images
if dx(2) >= 0
    xrange1 = [1+dx(2) size(im1, 2)];
    xrange2 = [1 size(im2,2)-dx(2)];
else
    xrange1 = [1-dx(2) size(im1,2)];
    xrange2 = [1 size(im2, 2)+dx(2)];
end

overlap_img1 = im1(yrange1(1):yrange1(2), xrange1(1):xrange1(2));
overlap_img2 = im2(yrange2(1):yrange2(2), xrange2(1):xrange2(2));

end