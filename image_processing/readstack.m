%%% read all images in a tiff stack into matlab
% updated 4/29/17
function imageStack = readstack(fname)

info = imfinfo(fname);
imageStack = [];
numberOfImages = length(info);
for k = 1:numberOfImages
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end 