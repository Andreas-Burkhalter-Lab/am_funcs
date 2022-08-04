%%% perform correlation analysis on intensity values across pixels in two
%%% images within a restricted part of the image

%%% 'file_area_mask' is a grayscale image containing nonzero values in the
%%% region to be analyzed and zeros in the region to be ignored
% all images must be the same dimensions
file_m2 = 'C:\Users\AM\Documents\Lab\sections\16041\16041_ctx_m2_s2c_x8_60s.tif'; 
file_egfp = 'C:\Users\AM\Documents\Lab\sections\16041\16041_ctx_egfp_s2c_x8_5s.tif'; 
file_area_mask = '16041_ctx_por-outline_s2c_x8.tif';
blur_radius = 15; % radius=15 for Rinaldo's patchiness analysis in Ji et al. 2015
show_image = 1; % display image of pixels within restricted area


% read images
m2 = double(imread(file_m2)); %
egfp = double(imread(file_egfp)); %% folder must be on matlab path
area_mask = logical(imread(file_area_mask));

% create blurred versions of images
blur_filter_2d = fspecial('disk',blur_radius); %%% create filter; arguments set blur radius
m2_blurred = imfilter(m2,blur_filter_2d,'replicate'); %% created blurred version
egfp_blurred = imfilter(egfp,blur_filter_2d,'replicate'); %% created blurred version

% obtain intensity values from the region of interest within each image
m2_roi = m2(area_mask);
egfp_roi = egfp(area_mask);
m2_blurred_roi = m2_blurred(area_mask);
egfp_blurred_roi = egfp_blurred(area_mask);

% perform correlation analysis
[roi_corrCoef_m2_egfp, roi_corrP_m2_egfp] = corrcoef(double(m2_roi),double(egfp_roi))
[roi_blurred_corrCoef_m2_egfp, roi_blurred_corrP_m2_egfp] = corrcoef(double(m2_blurred_roi),double(egfp_blurred_roi))
[wholeimage_corrCoef_m2_egfp, wholeimage_corrP_m2_egfp] = corrcoef(double(m2),double(egfp))
[wholeimage_blurred_corrCoef_m2_egfp, wholeimage_blurred_corrP_m2_egfp] = corrcoef(double(m2_blurred),double(egfp_blurred))

if show_image
    figure
    subplot(2,2,1)
    imagesc(m2.*double(area_mask))
    title('m2_roi')
    subplot(2,2,2)
    imagesc(egfp.*double(area_mask))
    title('egfp_roi')
    subplot(2,2,3)
    imagesc(m2_blurred.*double(area_mask))
    title('m2_blurred_roi')
    subplot(2,2,4)
    imagesc(egfp_blurred.*double(area_mask))
    title('egfp_blurred_roi')    
end
