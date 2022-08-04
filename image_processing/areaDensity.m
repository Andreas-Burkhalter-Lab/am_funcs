%%% get mean pixel values within delineated areas of series of input images  
%
%   outtable = areaDensity(intable, pars)
%
% input table or name of excel file with at least the following variables:
%   exposure_sec = exposure time; values will be converted to equivalent as if exposure time = 1s
%   proj_file = image with pixel values to be averaged
%   roi_rile = image with white areas indicating ROI, other areas black  
%   directory = folder in which all files for this row are found
%   baseline_file (optional).... see below
%   sub (optional) = subject; if included, normalization will use the max value of each subject, not of the whole table
%   
%   pars.area_proportion (default==1) takes only the top specified fraction of pixels from the ROI to analyze 
%       -(don't combine with intensity_vals_proportion)
%   pars.intensity_vals_proportion (default==1) takes only the top specified fraction of pixel intensity vals to analyze
%           intens_vals_prop_side (default='top')... 'bottom' or 'top'; which side of distribution to take fraction of pixels from
%   pars.do_baseline_subtraction (default==true) includes the step of baseline-subtracting intensity values
%   pars.separate_projfile_directory determines whether proj files are specified in a directory separate from roi and baseline files
%   pars.loadbw_ops specifies how black-white imges are loaded with loadbw.m
%
% three options for table variable 'baseline_file'
%       1. filename of .mat file with struct 'baseline' for normalizing to baseline (filename must have extension .mat):
%           baseline.baselineImage = filename of image to get baseline from
%           baseline.baseLineArea = filename of image in white indicating where to analyze baselineImage
%           baseline.exposure_sec = exposure time
%       2. filename of BW image indicating baseline area; use proj_file as intensity image  
%       3. no input; set pars.do_baseline_subtraction to false to skip baseline subtracting
%
% only file paths listed in the 'directory' variable will be considered; paths in other table variables will be ignored
%
%%%% last upated 2021/11/1 on thermaltake

function outtable = areaDensity(intable, pars)

%% set options
show_plot = 0; 
    sectionthickness = 40; % microns.... only affects plotting, not data analysis

pars = vardefault('pars', struct);
pars.blur_image = field_default(pars, 'blur_image', 0); % choose whether or not to blur before anlayzing pixels
    pars.blur_using_only_roi = field_default(pars,'blur_using_only_roi',1); % if true default), do not use non-roi pixels for blurring 
    pars.diskblurradius_um = field_default(pars,'diskblurradius_um',29); % if blurring, this will be radius of filter
pars.save_analyzed_image = field_default(pars,'save_analyzed_image',0); % if true, save the analyzed image in the output table (may be useful for looking at blurred images)
pars.area_proportion = field_default(pars, 'area_proportion', 1); % if not specified, analyze the whole ROI (proportion 1/1)
pars.intensity_vals_proportion = field_default(pars, 'intensity_vals_proportion', 1); % if not specified, analyze all intensity vals (proportion 1/1)
        pars.intens_vals_prop_side = field_default(pars,'intens_vals_prop_side','top'); % default to taking brightest pixels
pars.do_baseline_subtraction = field_default(pars, 'do_baseline_subtraction', true); % if not specified, perform baseline subtraction
pars.separate_projfile_directory = field_default(pars, 'separate_projfile_directory', false); % if true, read full filepath for proj files, don't use directory var
pars.output_analyzed_pixels = field_default(pars,'output_analyzed_pixels',false); % if true, output all analyzed pixels for each case in outtable... uses a lot of memory
pars.loadbw_ops = field_default(pars,'loadbw_ops',struct);

%%
if ischar(intable) % intable assumed to be excel filename
    filetable = readtable(intable);
else % intable assumed to be table variable
    filetable = intable;
end
nimages = height(filetable);
outtable = filetable;
outtable.roitotalpix = NaN(nimages,1);
outtable.intensPerPix = NaN(nimages,1);

%%%% if exposure time isn't specified, default to 1
if any(strcmp('exposure_sec', filetable.Properties.VariableNames)) % if not specified
    outtable.exposure_sec = ones(nimages,1); % default to 1
    fprintf('     Exposure times not specified; defaulting to 1\n')
end

%% analyze each image
wbar = waitbar(0,'Analyzing intensity in specified area...');
for ifile = 1:nimages
    %%% load proj file and roi file from the directories listed
    if pars.separate_projfile_directory
        projFileToLoad = filetable.proj_file{ifile}; % don't assume proj file is in same directory as roi and baseline files
    elseif ~pars.separate_projfile_directory
        projFileToLoad = [filetable.directory{ifile}, filesep, getfname(filetable.proj_file{ifile}, 'extension')];
    end
    projimg = imread(projFileToLoad);
    projimg = projimg(:,:,1); % in case image is rgb, take only first channel
    roiFileToLoad = [filetable.directory{ifile}, filesep, getfname(filetable.roi_file{ifile}, 'extension')];
    roiimg = loadbw(roiFileToLoad, pars.loadbw_ops);
    
    if any(size(projimg) ~= size(roiimg))
        error([outtable.proj_file{ifile} ' and ' outtable.roi_rile{ifile} ' are different sizes.'])
    end

    projimgroi = projimg; % projimgroi will be projimg but with pixels outside of the ROI turned into NaN
    projimgroi = double(projimgroi);  % must be double to have nans
    projimgroi(~roiimg) = NaN; % set non-roi pixels to nan
    
    % blurring 
    if pars.blur_image
        diskblurradius_pix = pars.diskblurradius_um * pixPerUm(filetable.zoom(ifile),filetable.scope{ifile});
        H = fspecial('disk',diskblurradius_pix); %%% create filter; arguments set blur radius

        if pars.blur_using_only_roi
            imroiblurred = nanconv(projimgroi,H,'nanout','edge');
        else
            imblurred = imfilter(im,H,'replicate'); %% 
            imroiblurred = imblurred;
            imroiblurred(~roiimg) = NaN;
        end
        projimgroi = imroiblurred; 
    end
    
    % get ROI intensity vals
    roi_vals = projimgroi(roiimg);
    if pars.area_proportion ~= 1 % if we are using only a fraction of the pixels from ROI
        roi_vals_sorted = flipud(sort(roi_vals)); % highest values listed first
        nvalues_to_use = round(pars.area_proportion * length(roi_vals)); 
        roi_vals = roi_vals_sorted(1:nvalues_to_use); % cut down ROI values to those in the top fraction
    end
    if pars.intensity_vals_proportion ~=1 % if we are using only a fraction of the intensity vals
        if pars.area_proportion ~= 1 
            error('Do not use pars.area_proportion and pars.intensity_vals_proportion together')
        end
        if strcmp(pars.intens_vals_prop_side, 'top')
            intensity_vals_sorted = sort(unique(roi_vals),'descend');  % highest values listed first
            n_unique_intensity_vals_to_use = round(pars.intensity_vals_proportion * length(intensity_vals_sorted));
            outtable.intensity_threshold_raw(ifile) = intensity_vals_sorted(n_unique_intensity_vals_to_use); % discard pixel intensity values below this threshold
            roi_vals = roi_vals(roi_vals >= outtable.intensity_threshold_raw(ifile)); % apply threshold
        elseif strcmp(pars.intens_vals_prop_side, 'bottom')
            intensity_vals_sorted = sort(unique(roi_vals),'ascend');  % lowest values listed first
            n_unique_intensity_vals_to_use = round(pars.intensity_vals_proportion * length(intensity_vals_sorted));
            outtable.intensity_threshold_raw(ifile) = intensity_vals_sorted(n_unique_intensity_vals_to_use); % discard pixel intensity values above this threshold
            roi_vals = roi_vals(roi_vals <= outtable.intensity_threshold_raw(ifile)); % apply threshold
        else
            error('ops.intens_vals_prop_side must be set to either ''top'' or ''bottom''')
        end

    end
    if pars.output_analyzed_pixels % output all analyzed pixels for each case in outtable
        outtable.analyzed_pix{ifile} = roi_vals;
    end
    raw_roi_mean = nanmean(roi_vals);
    
    % baselining
    if pars.do_baseline_subtraction  % do baseline subtraction
        baseFileToLoad = [filetable.directory{ifile}, filesep, getfname(filetable.baseline_file{ifile}, 'extension')];
        if strcmp(baseFileToLoad(end-3:end), '.mat')   % if baseline input is a .mat file with baseline struct
            baseline = load(baseFileToLoad,'baseline'); baseline = baseline.baseline; % load baseline struct for use with this image
            baseimg = imread(baseline.baselineImage);
            baseimg = baseimg(:,:,1); % in case image is rgb, take only first channel
            baseAreaImg = loadbw(baseline.baselineAreaImage, pars.loadbw_ops);
            baseVals = baseimg(baseAreaImg);
            outtable.baseline(ifile) = nanmean(baseVals) / baseline.exposure_sec;
            outtable.intensPerPixRaw(ifile) = raw_roi_mean / outtable.exposure_sec(ifile); % record exposure-normed raw val before subtracting baseline
            outtable.intensPerPix(ifile) = outtable.intensPerPixRaw(ifile) - outtable.baseline(ifile);
        elseif ~strcmp(baseFileToLoad(end-3:end), '.mat')    % if baseline input is a BW image file
            baseAreaImg = loadbw(baseFileToLoad, pars.loadbw_ops);
            projimgbase = projimg; % projimgbase will be projimg but with pixels outside of the baseline image turned into NaN
            projimgbase = double(projimgbase);  % must be double to have nans
            projimgbase(~baseAreaImg) = NaN; % set non-baseimg pixels to nan
            
            if pars.blur_image % blur basline area image
                if pars.blur_using_only_roi
                    projimgbase = nanconv(projimgbase,H,'nanout','edge');
                else
                    imblurred = imfilter(im,H,'replicate'); %% 
                    projimgbase = imblurred;
                    projimgbase(~baseAreaImg) = NaN;
                end
            end

            baseVals = projimgbase(baseAreaImg);
            outtable.intensPerPixRaw(ifile) = raw_roi_mean / outtable.exposure_sec(ifile); % record exposure-normed raw val before subtracting baseline
            raw_roi_mean = raw_roi_mean - nanmean(baseVals); % subtract baseline
            outtable.baseline(ifile) = nanmean(baseVals) / outtable.exposure_sec(ifile); 
            outtable.intensPerPix(ifile) = raw_roi_mean / outtable.exposure_sec(ifile); % normalize to simulated exposure time of 1s
        end
    elseif ~pars.do_baseline_subtraction % don't do baseline subtraction
        outtable.intensPerPix(ifile) = raw_roi_mean / outtable.exposure_sec(ifile); % normalize to simulated exposure time of 1s
    end

    outtable.roitotalpix(ifile) = nnz(roiimg);
    outtable.roitotalpix_used(ifile) = length(roi_vals); % how many pixels of the ROI were used in the analysis
    if pars.save_analyzed_image
        projimg = double(projimg);  % must be double to have nans
        projimg(~roiimg) = NaN; 
        outtable.analyzed_image{ifile} = projimgroi; % output the entire analyzed image, after 
    end
    if isvalid(wbar)
        waitbar(ifile/nimages,wbar)
    end
end
if isvalid(wbar)
    close(wbar)
end    

%% normalization of intensity values to maximum value within each case
if any(strcmp('sub',filetable.Properties.VariableNames))  %%% if a variable for subject is included in the table, norm using only this subject
    subs = unique(outtable.sub);
    for isub = 1:length(subs)
        thissub = subs(isub);
        these_rows = outtable.sub==thissub;
        outtable.normedIntens(these_rows) = outtable.intensPerPix(these_rows) / max(outtable.intensPerPix(these_rows));
    end
else % norm using the entire table
     outtable.normedIntens = outtable.intensPerPix / max(outtable.intensPerPix);
end
    
    
    
    
    
    
    %% plotting
if show_plot && any(strcmp(outtable.Properties.VariableNames,'sec'))
    figure
    secloc = sectionthickness * outtable.sec;  %%%%% - sectionthickness/2;
    plot(secloc,outtable.normedIntens)
    xlabel('Depth (microns)')
    ylabel('Projection Intensity Per Pixel')
%     set(gca,'ylim',[0 max(get(gca,'ylim'))])
end



