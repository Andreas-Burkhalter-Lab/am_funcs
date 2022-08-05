function [images, masks] = jpm_region_mask(entry)
% jpm_region_mask(entry) returns a cell of masked images from an xdb entry
%     Syntax: images = jpm_region_mask(entry) 
%     Purpose: This function takes an entry with the 'hist_analysis.regions'
%             subfield and returns a cell variable 'images' with the same number of
%             elements as data_locations.
%             
%             Each cell of images will, in turn, contain a number of masked
%             images based on the number of defined regions in
%             hist_analysis.
%             
%             Thus, if you have 3 data_locations (tif images) in an entry, each of which 
%                  contains 2 defined regions, your cell layout will be:
%                  [ [1x2 cell] [1x2 cell] [1x2 cell] ] and each element in
%                  this variable will contain an image with the same base
%                  size as the image in each data_locations (.tif) file
%
%     See also IMREAD, IMSHOW, POLY2MASK, CELLCOUNT_FLAT

%  Copyright 2007 Julian P. Meeks (Timothy Holy Laboratory)
%  
%  Revision History: 10/10/2007: wrote it (JPM)

% Step 1: count the number of images (from entry.data_locations)
if ~isfield(entry, 'data_locations')
    images = []; return;
elseif isempty(entry.data_locations)
    images = []; return;
end
n_images = size(entry.data_locations);
images = cell(n_images);                 % set up <images> variable
masks = cell(n_images);

% Step 2: count the number of regions for each image
if ~isfield(entry, 'hist_analysis')
    images = []; return;
elseif ~isfield(entry.hist_analysis, 'regions')
    images = []; return;
end
for idx = 1:n_images(2)
    n_regions(idx) = size(entry.hist_analysis.regions{idx},2);
end

% Step 3: transform region points into x,y coordinates
region_xpts = cell(1); region_xpts(1) = [];
region_ypts = cell(1); region_ypts(1) = [];
for idx_images = 1:n_images(2)
    for idx_regions = 1:n_regions(idx_images)
        xpts = [];
        ypts = [];
        for idx_region_pts = 1:size(entry.hist_analysis.regions{idx_images}{idx_regions},2)
            xpts(idx_region_pts) = entry.hist_analysis.regions{idx_images}{idx_regions}{idx_region_pts}(1);
            ypts(idx_region_pts) = entry.hist_analysis.regions{idx_images}{idx_regions}{idx_region_pts}(2);
        end
        these_xpts(idx_regions) = {xpts};
        these_ypts(idx_regions) = {ypts};
    end
    region_xpts(idx_images) = {these_xpts};
    region_ypts(idx_images) = {these_ypts};
    these_xpts(:,:) = [];
    these_ypts(:,:) = [];
end

% Step 4:create the masked images
for idx_images = 1:n_images(2)
    base_image = double(imread(entry.data_locations{idx_images}));
    for idx_regions = 1:n_regions(idx_images)
            these_masks(idx_regions) = {poly2mask(region_xpts{idx_images}{idx_regions},...
                                               region_ypts{idx_images}{idx_regions},...
                                               size(base_image,1),...
                                               size(base_image,2))};
            these_images(idx_regions) = {base_image.*these_masks{idx_regions}};
    end
    images(idx_images) = {these_images};
    masks(idx_images) = {these_masks};
end


end