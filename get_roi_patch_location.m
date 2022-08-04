%%%% get_roi_patch_location: determine whether rois are in patches or not
% tuning_struct_in = results struct from get_tuning_curves
% patchdata = output struct  from analyzePatches
% registration_file = file created when warping in vivo time-averaged image
%       to match m2 image (or gcamp image to which m2 image is registered)
%%%
%%% tuning_struct_out = get_roi_patch_location(tuning_struct_in,patchdata_file,registration_file)
%%%
%%% updated 2019/1/28 on thermaltake

function tuning_struct_out = get_roi_patch_location(tuning_struct_in, patchdata_file, registration_file)

tuning_struct_out = tuning_struct_in; 
reg_struct = load(registration_file); 
tuning = tuning_struct_in.tuning;
nrois = height(tuning);
tuning.inpatch = false(nrois,1);
tuning.patch_center_dist = NaN(nrois,1);
tuning.patch_dist = NaN(nrois,1);
tuning.nearest_patch_center_yx = NaN(nrois,2);
tuning.Properties.VariableUnits(strcmp(tuning.Properties.VariableNames,'patch_center_dist')) = {'um'};
tuning.Properties.VariableNames(strcmp(tuning.Properties.VariableNames,'centeryx')) = {'centeryx_prereg'}; % clarify that these points are not registered to m2 image yet
tuning.Properties.VariableNames(strcmp(tuning.Properties.VariableNames,'cellimage')) = {'cellimage_prereg'}; % clarify that these points are not registered to m2 image yet
tuning.centeryx_reg = NaN(nrois,2); % centers of rois projected into m2 image
tuning.roi_image_reg = cell(nrois,1); 
patchdata = load(patchdata_file); patchdata = patchdata.patchdata;
patchimage = patchdata.patchimage;
if isfield(patchdata,'patch_centers_image'); [patch_centers_yx(:,1), patch_centers_yx(:,2)] = find(patchdata.patch_centers_image); end
[patch_pix_yx(:,1), patch_pix_yx(:,2)] = find(patchimage);

% images size checks
if any(size(tuning.roi_image_prereg{1}) ~= size(reg_struct.movingimage))
    error('moving image dimension mismatch')
end
refimage = reg_struct.refimage(:,:,1);
if any(size(patchdata.im) ~= size(refimage))
    error('reference image dimension mismatch')
end

%%% project rois onto m2 image, get new registered points
roicenters_reg = NaN(size(refimage));
for iroi = 1:nrois
%     iroi
    thisregmask = imwarp(full(tuning.roi_image_prereg{iroi}),reg_struct.geotransform,'OutputView',reg_struct.worldref);
    tuning.roi_image_reg{iroi} = sparse(thisregmask); 
    if any(thisregmask(:)) % if this roi is still in frame after warping
        [rows, cols] = ndgrid(1:size(thisregmask,1), 1:size(thisregmask,2));
        tuning.centeryx_reg(iroi,1) = round(sum(sum(rows(thisregmask))) ./ sum(sum(thisregmask))); % y center from top
        tuning.centeryx_reg(iroi,2) = round(sum(sum(cols(thisregmask))) ./ sum(sum(thisregmask))); % x center from left
        roicenters_reg(tuning.centeryx_reg(iroi,1), tuning.centeryx_reg(iroi,2)) = iroi; % add point to image of reg'd centers
    else % if the roi got warped out of frame
        warning(['Roi ' num2str(iroi) ' got registered to a location out of the frame and will not be analyzed for patch-relative location.'])
    end
end

%%%% find roi locs relative to patches and quantiles
for iroi = 1:nrois
    if ~isnan(tuning.centeryx_reg(iroi,1))
        tuning.inpatch(iroi) = patchimage(tuning.centeryx_reg(iroi,1),tuning.centeryx_reg(iroi,2));
        if isfield(patchdata,'patch_centers_image')
            [nearestcenter, tuning.patch_center_dist(iroi)] = dsearchn(patch_centers_yx,tuning.centeryx_reg(iroi,:));
            tuning.nearest_patch_center_yx(iroi,:) = patch_centers_yx(nearestcenter,:);
        elseif any(ismember(tuning.Properties.VariableNames,'patch_center_dist'))
            tuning.patch_center_dist = [];
            tuning.nearest_patch_center_yx = [];
        end
        [~,tuning.patch_dist(iroi)] = dsearchn(patch_pix_yx,tuning.centeryx_reg(iroi,:)); % distance from any patch pix; ==0 if in patch
        tuning.quantile(iroi) = patchdata.quant_levels_img(tuning.centeryx_reg(iroi,1),tuning.centeryx_reg(iroi,2));
    end
end

tuning_struct_out.tuning = tuning; 
tuning_struct_out.registration_file = registration_file; 
tuning_struct_out.reg_struct = reg_struct; 
tuning_struct_out.patchdata_file = patchdata_file;
tuning_struct_out.patchdata = patchdata; 