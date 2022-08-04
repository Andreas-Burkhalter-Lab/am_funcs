function [imageout imfile] = loadDensityImage(imagein)

if ischar(imagein)
    imfile = imagein;
    imagein = imread(imagein);
else
    imfile = '';
end

imagein = imagein(:,:,1);
imageout = imagein;