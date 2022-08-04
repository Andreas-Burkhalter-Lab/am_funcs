%%%% record the xy location and center of an roi into the tuning table and rotate if appropriate               
% updated 2018/12/04 on thermaltake

% rotate han lab image so that anterior is upward; anterior should already be upward on 4th floor scope recordings
if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
    roiimage_pre_rotate = sparse(squeeze(logical(dff_data.masks(iroi,:,:))));
    tuning.roi_image_pre_rotate{iroi} = roiimage_pre_rotate;
    roi_image_prereg = rot90(roiimage_pre_rotate,-1); %%% cell mask after rotating to make anterior point upward
    [rows, cols] = ndgrid(1:size(roiimage_pre_rotate,1), 1:size(roiimage_pre_rotate,2));
    tuning.centeryx_pre_rotate(iroi,1) = round(sum(sum(rows(roiimage_pre_rotate))) ./ sum(sum(roiimage_pre_rotate))); % y center from top
    tuning.centeryx_pre_rotate(iroi,2) = round(sum(sum(cols(roiimage_pre_rotate))) ./ sum(sum(roiimage_pre_rotate))); % x center from left
elseif strcmp(stimpars.computer,'ANDREWLAB-PC')  %% if recorded on 4th floor scope, don't rotate
    roi_image_prereg = sparse(squeeze(logical(dff_data.masks(iroi,:,:))));
end
tuning.roi_image_prereg{iroi} = roi_image_prereg;

% get center of mass of cell image - use this location as an identifier for the cell
[rows, cols] = ndgrid(1:size(roi_image_prereg,1), 1:size(roi_image_prereg,2));
tuning.centeryx(iroi,1) = round(sum(sum(rows(roi_image_prereg))) ./ sum(sum(roi_image_prereg))); % y center from top
tuning.centeryx(iroi,2) = round(sum(sum(cols(roi_image_prereg))) ./ sum(sum(roi_image_prereg))); % x center from left