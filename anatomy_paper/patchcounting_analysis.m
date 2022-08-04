%%%% for each of v1, por, lm, li, find m2 patches in example images and compare spatial distribution and shape features
% [cases, patchtable] = patchcounting_analysis(patchcounting_filelist,excel_rows)
%
%%%%%%% input should be excel spreadsheet with columns labeled: area, sub, m2_file, roi_file, roi_file_for_area_orient, zoom, scope
%
% updated 20/2/18

function [cases, patchtable] = patchcounting_analysis(patchcounting_filelist,excel_rows,findpatches_pars)

plot_patch_shape = 0; 

% magnification factor data from Garrett et al. 2014, fig 7c
magfactor_data = table({'v1';'lm';'li';'por'},[9.4;7.8;4.89;6.7]*1e-4,'VariableNames',{'area','sqmm_per_sqdeg'}); 

%% import file list
cases = readtable(patchcounting_filelist); 
excel_rows = vardefault('excel_rows',2:1+height(cases)); % default to using all cases
cases = cases(excel_rows-1,:); % only include the specified case rows
nfiles = height(cases);

%% analyze patch spatial distrubution for each case
findpatches_pars = vardefault('findpatches_pars',struct); % make the pars struct if it wasn't specified
barhandle = waitbar(0,['Analyzing patches...']);
for ifile = 1:nfiles
    cases.patchdata{ifile} = findpatches(cases.m2_file{ifile},cases.roi_file{ifile},cases.zoom(ifile),cases.scope{ifile},findpatches_pars);
    cases.density{ifile} = patchDensity(cases.patchdata{ifile}.patchimage,cases.patchdata{ifile}.roi,cases.zoom(ifile),cases.scope{ifile});
    cases.meandistum(ifile) = cases.density{ifile}.meandistum;
    cases.patches_per_sqmm(ifile) = cases.density{ifile}.patches_per_sqmm; %% areal patch density
     if any(strcmp(magfactor_data.area,cases.area{ifile})) % if this area is listed in magfactor table
         this_area_sqmm_per_sqdeg = magfactor_data.sqmm_per_sqdeg(strcmp(magfactor_data.area,cases.area{ifile})); % magfactor for this area
         cases.patches_per_sqdeg(ifile) =  this_area_sqmm_per_sqdeg * cases.patches_per_sqmm(ifile); %% compute retinotopic patch density from areal density and magnification factor
     end
    cases.density{ifile} = rmfield(cases.density{ifile},{'patchimage','roi','patchesimage_roi'});
    roi_for_area_orient = loadbw(cases.roi_file_for_area_orient{ifile}); % load image for determing major axis orientation of area
    cases.roi_orient(ifile) = cell2mat(struct2cell(regionprops(roi_for_area_orient,'Orientation'))); % get angle of major axis of roi
    deg_to_rotate = 90 - cases.roi_orient(ifile); 
    % get  matrix for for rotating patches so that y axis represents major axis of area, x axis represents minor axis of area (use negative area orientation)
    rotation_mat = [cosd(deg_to_rotate) -sind(deg_to_rotate); sind(deg_to_rotate) cosd(deg_to_rotate)]; 
    % get length of area along its own major and minor axes
    [y, x] = find(roi_for_area_orient); 
    arearotated = round(rotation_mat * [y x]')';
    area_length_um_major = range(arearotated(:,1)) * umPerPix(cases.zoom(ifile), cases.scope{ifile}); % length along y axis after rotation
    area_length_um_minor = range(arearotated(:,2)) * umPerPix(cases.zoom(ifile), cases.scope{ifile}); % length along x axis after rotation
    cases.area_major_minor_ratio(ifile) = area_length_um_major / area_length_um_minor; 
    cases.patchstats{ifile} = table;
    npatches = length(cases.patchdata{ifile}.patch_bnds); % ? causes errors with certain findpatches parameters? quantiles? 
    for ipatch = 1:npatches
        patch_bnds_rotated_yx = [rotation_mat * cases.patchdata{ifile}.patch_bnds{ipatch}']'; % rotate this patch so that 0 deg = aligned with area major axis
        cases.patchstats{ifile}.length_um_major(ipatch) = range(patch_bnds_rotated_yx(:,1)) * umPerPix(cases.zoom(ifile), cases.scope{ifile}); % length along y axis after rotation
        cases.patchstats{ifile}.length_um_minor(ipatch) = range(patch_bnds_rotated_yx(:,2)) * umPerPix(cases.zoom(ifile), cases.scope{ifile}); % length along x axis after rotation
        cases.patchstats{ifile}.major_minor_ratio(ipatch) = cases.patchstats{ifile}.length_um_major(ipatch) / cases.patchstats{ifile}.length_um_minor(ipatch); % major:minor axis length ratios
    end
    cases.mean_patch_major_minor_ratio(ifile) = mean(cases.patchstats{ifile}.major_minor_ratio);
    cases.med_patch_major_minor_ratio(ifile) = median(cases.patchstats{ifile}.major_minor_ratio);
    try waitbar(ifile/nfiles,barhandle);end
end    
try close(barhandle); end
cases = movevars(cases, 'area_major_minor_ratio', 'Before', 'scope');
cases = movevars(cases, 'mean_patch_major_minor_ratio', 'Before', 'area_major_minor_ratio');
cases = movevars(cases, 'patches_per_sqmm', 'Before', 'area_major_minor_ratio');

%% make table of all patches to compare across areas
patchtable = table;
for areaname = unique(cases.area)'
    for ii = find(strcmp(cases.area,areaname))'
        patchtable = [patchtable; [table(repmat(areaname,height(cases.patchstats{ii}),1),'VariableNames',{'area'}), cases.patchstats{ii} ]];
    end
end

%% plotting
if plot_patch_shape
%     [anovap, anovatab, anovastats] = anova1(patchtable.major_minor_ratio, patchtable.area);
    [anovap, anovatab, anovastats] = anova1(cases.mean_patch_major_minor_ratio, cases.area);
end
