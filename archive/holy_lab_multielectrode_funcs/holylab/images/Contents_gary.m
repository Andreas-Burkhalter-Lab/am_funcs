ai_read % This script is used for reading the .ai files containing the analogsignals stored by Imagine
analyze_example.m %  Example analysis script
brick_neighbors.m % Find the bricks that share vertices (useful in breaking up multidimenional clusters).
calc_stack_roi_inten.m % Takes roi definitions within stack space, created by stack_roi,
%and calculates the mean pixel intensity within each of those regions over time. 
calc_stack_roi_inten_new.m % Takes roi definitions within stack space, created by stack_roi,
% and calculates the mean pixel intensity within each of those regions over time.
calculate_region_area.m % CALCULATE_REGION_AREA calculates the area of a given region from
% hist_analysis field in .xdb database
cameradata.m % a structure for specifying camera performance where I (intensity) is related to I = gain * n + eta + bias
cellcount_flat.m % GUI interface for cell counting within regions of flattened images
clust_roi_explore.m % connect responses after clustering to raw images back to "raw" images
colorize_dfof.m % create colorized images to encode DeltaF/F
colorize_pixels.m % generate a color image from grayscale, using color channel masking
color_marked_pixels_red.m % display image in grayscale except chosen pixels are red
convert_stack_to_avi.m % converts a basic stack to an avi movie

==> Deconvolution of images which have a spatial sampling frquency greater than the light sheet thickness
deconv_exponential.m % estimate underlying signal from filtered, noisy version
deconv_filter_create.m % create a Weiner deconvolution filter for one-dimensional data
deconv_gaussian.m % create a deconvolution filter for gaussian blur
deconvstack_design.m % viewing tool for 4D image stacks
deconvstack_execute.m % no help


find_bad_frames.m % Compute rigid registration and mismatch on a per-frame basis for whole movie 
flattentifstack.m % opens a gui-led process to flatten 16-bit tif stacks
frame_error.m % compute mismatch on a per-frame basis for whole movie

gaussian_filter.m % create multidimensional gaussian filters
gaussians2im.m % calculate a model image, given gaussian parameters
im2gaussians.m % find gaussian-shaped objects in an image
image_expand.m % image_expand takes a small image and expands the exact image to a larger size (for any reason you might like)
image_shift.m % create a translated image
image_snip_qinterp.m % cut out rectangular regions of an image, with subpixel interpolation
imagine2med.m % convert imagine files to medical image file formats: Mayo Analyze(.img), NIfTI(.nii), and DICOM(.dcm) files.
imagine2multitiff.m % takes .imagine binary data and saves it down in multi-tiff format. 
