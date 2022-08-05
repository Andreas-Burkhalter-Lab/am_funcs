function reorient_roidef_file(roidef_file)

% this function will take in an roi_def file and add a new field which will
% consist of the true (tissue) coordinates of the position in um

% Example usage:
%  reorient_roidef_file('vno_gcamp2_2010_01_18_2_reg_pixel_corrected.roidef')
% 
% This will add a new field called posInUm_tissue
% This is the position in tissue coordinates. This assumes imaging at 45
% degrees
%
% Diwakar Turaga 2010

load(roidef_file,'-mat');

theta = 45*pi/180; 
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for ROI_indx = 1:size(roi_defs,2) 
    roi_pos  = roi_defs(ROI_indx).posInUm;
    roi_pos = roi_pos';
    roi_pos(2:3,:) = R*roi_pos(2:3,:);
    roi_pos = roi_pos';
    roi_defs(ROI_indx).posInUm_tissue = roi_pos;
    
end

save(roidef_file, 'pixelPerUm', 'roi_defs','header');
