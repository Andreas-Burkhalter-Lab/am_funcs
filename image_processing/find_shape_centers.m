%%% get centers of white objects from a bw image within a specified roi
% [shapeCentersInds_raw, shapeCentersInds_roiRelative, shapeCentersImage] = find_shape_centers(bw_shapes_image, roi_image)
%%%% updated 18/6/4 on thermaltake
  
function [shapeCentersInds_raw, shapeCentersInds_roiRelative, shapeCentersImage] = find_shape_centers(bw_shapes_image, roi_image)
  
[bounds, labelmat, nshapes_total, adjc] = bwboundaries(bw_shapes_image);
shapeCentersYX_raw = NaN(nshapes_total,2);
for j = 1:nshapes_total
    [y, x] = find(labelmat == j);
    shapeCentersYX_raw(j,1) = round(mean(y));
    shapeCentersYX_raw(j,2) = round(mean(x));
end
shapeCentersInds_raw = sub2ind(size(bw_shapes_image),shapeCentersYX_raw(:,1),shapeCentersYX_raw(:,2));
%%% remove shapes outside of the roi
if exist('roi_image','var')
    roiinds = find(roi_image);
    shapeCentersInds = intersect(shapeCentersInds_raw, roiinds);
    shapeCentersImage = false(size(bw_shapes_image));
    shapeCentersImage(shapeCentersInds) = true; 
    shapeCentersInds_roiRelative = shapeCentersImage(roi_image); % get indices of shape centers WITHIN the roi, as a logical vector
end