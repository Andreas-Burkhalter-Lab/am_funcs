%%%% FINDPATCHES: automatically find m2 patches in an image based on intensity peaks
% % patchdata = findpatches(imfile,roifile,zoom,scope,pars)
% %  pars includes: minAreaPatch_squm
            %       maxAreaPatch_squm
            %       diskblurradius_um
            %       threshMode
            %       threshFraction
%%% updated 2020/10/1 on thermaltake
 
 function patchdata = findpatches(imfile,roifile,zoom,scope,pars)
pars = vardefault('pars',struct);

%% get pars, set defaults
% min area in square microns a blob must contain to be considered a patch
pars.minAreaPatch_squm = field_default(pars,'minAreaPatch_squm',0);  
%     minAreaPatch_squm = 236; % 236 seems best for subject 09060 for replicating patch density values from ji et al 2015 (patches top 2/6 quantiles)

% break up patches that are larger than this size in square microns
pars.maxAreaPatch_squm = field_default(pars,'maxAreaPatch_squm',inf); 
% pars.maxAreaPatch_squm = field_default(pars,'maxAreaPatch_squm',7000);
% pars.maxAreaPatch_squm = field_default(pars,'maxAreaPatch_squm',2300);% 2300 seems best for subject 09060 for replicating patch density values from ji et al 2015

pars.blur_using_only_roi = field_default(pars,'blur_using_only_roi',1); % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
pars.diskblurradius_um = field_default(pars,'diskblurradius_um',29); % radius in microns of circular averaging filter

% choose method for picking the intensity threshold
pars.threshMode = field_default(pars,'threshMode','intensityQuantiles'); % divide pixels in discreet quantiles based on intensity
    pars.nQuantiles = field_default(pars,'nQuantiles',6); % number of intensity levels to divide image into (if using 'intensityQuantiles')
    pars.patchQuantiles = field_default(pars,'patchQuantiles',[4 5 6]); % % quantiles selected to be counted as patches(if using 'intensityQuantiles'); higher numbers = brighter pixels
% pars.threshMode = field_default(pars,'threshMode','fractionOfPopulation');
%     pars.fracPixSuperthresh = field_default(pars,'fracPixSuperthresh',0.5); % fraction of pixels marked as above threshold (if using 'fractionOfPopulation')
% pars.threshMode = field_default(pars,'threshMode','fractionOfMaxMinusMin');
%     pars.threshFraction = field_default(pars,'threshFraction',0.8); %%% thresh needed to put a pixel within a patch; value = fraction of intensity distance from min val within roi to max val within roi
% pars.threshMode = field_default(pars,'threshMode','fractionOfMax');
%     pars.threshFraction = field_default(pars,'threshFraction',0.878);

% if include_interior_nonroi_in_roi==true, include non-roi pixels surrounded by roi within the highest adjacent quantile
pars.include_interior_nonroi_in_roi = field_default(pars,'include_interior_nonroi_in_roi',1); 
pars.show_plots = field_default(pars,'show_plots',1);
pars.visualize_shrinking = field_default(pars,'visualize_shrinking',0); % plot the results of changing the threshold
pars.raisethresh_increment = field_default(pars,'raisethresh_increment',0.001); % how much to increase threshold each iteration when finding thresh; may get stuck if below 0.3 on x16 zoom
pars.save_patchdata = field_default(pars,'save_patchdata',1); % save patchdata struct and quants image after doing analysis
pars.quantborders_ops = field_default(pars,'quantborders_ops',struct); % ops for save_quant_borders_image
pars.loadbw_ops = field_default(pars,'loadbw_ops',struct); % ops for loadbw.m used for load roi file

scope = vardefault('scope','epimicro');




%%
if strcmp(pars.threshMode, 'intensityQuantiles') && pars.minAreaPatch_squm > 0 || pars.maxAreaPatch_squm < inf
        error('do not use patch size limits with intensityQuantiles')
end

minAreaPatch_pix = pars.minAreaPatch_squm * (pixPerUm(zoom,scope))^2;
maxAreaPatch_pix = pars.maxAreaPatch_squm * (pixPerUm(zoom,scope))^2;
diskblurradius_pix = pars.diskblurradius_um * pixPerUm(zoom,scope);
patchdata.minAreaPatch_pix = minAreaPatch_pix; 
patchdata.maxAreaPatch_pix = maxAreaPatch_pix;
patchdata.diskblurradius_pix = diskblurradius_pix; 

%% load images
 if ischar(imfile)
    im = imread(imfile);
else
    im = imfile;
    imfile = ' ';
end
im = im(:,:,1); % in case image is rgb, take only first channel
 
if exist('roifile','var') && ~isempty(roifile)
    if ischar(roifile)
        save_prepend = strrep(getfname(roifile), 'roi', '');
    end
    [roi roifile] = loadbw(roifile,pars.loadbw_ops);
    patchdata.roifile = roifile;
    roi = roi>0; 
elseif ~exist('roifile','var') || isempty(roifile) % analyze all pixels
    patchdata.roifile = [];
    roi = true(size(im));
end
save_prepend = vardefault('save_prepend','_'); 

if any(~(size(roi) == size(im)))
    error('image and roi file dimension mismatch')
end

findpatches_patchbounds(); %%%% find threshold intensity for patches

%% apply patch size limits
% break up large patches
if pars.maxAreaPatch_squm < inf
    blobs_labels = unique(labelmat_allsizes);
    blobs_labels(blobs_labels==0) = [];
    bnds_allsizes_presplit = bnds_allsizes;
    labelmat_allsizes_presplit = labelmat_allsizes;
    for i = 1:length(blobs_labels)
        blobsizes(i) = length(find(labelmat_allsizes==blobs_labels(i)));
    end
    maxblobsizefound = max(blobsizes);

    barhandle = waitbar(0,['Finding patches in ' imfile '...']);


    while maxblobsizefound > maxAreaPatch_pix
        i = find(blobsizes>maxAreaPatch_pix,1); % find blob to split
        thisbloblabel = blobs_labels(i);
        thisblobimage = imroiblurred;
        thisblobimage(labelmat_allsizes ~= thisbloblabel) = 0;
        npixnew = length(find(labelmat_allsizes==thisbloblabel)); 
        newthresh = thresh;
        nblobshere = 1;
        if length(blobs_labels) > 5
            waitbar(thisbloblabel/length(blobs_labels),barhandle);
        end
        while npixnew > maxAreaPatch_pix && nblobshere < 2
            newthresh = newthresh + threshinc;
            superthreshhere = thisblobimage>=newthresh;
            [newbnds, newlabelmat, nblobshere] = bwboundaries(superthreshhere,'noholes');
            npixnew = length(find(superthreshhere));
            if pars.visualize_shrinking
                hold off; imagesc(im);       hold all;        scatter(newbnds{1}(:,2),newbnds{1}(:,1),'r','.');     colorbar; pause(0.1);
            end
        end
        newlabelmat(newlabelmat>1) = newlabelmat(newlabelmat>1) + length(bnds_allsizes); % relabel new blobs
        %%% add new blobs to image and list
        labelmat_allsizes(labelmat_allsizes==thisbloblabel) = 0;
        labelmat_allsizes(newlabelmat==1) = thisbloblabel;    
        bnds_allsizes{i} = newbnds{1};
        if nblobshere > 1
            bnds_allsizes = [bnds_allsizes; newbnds(2:end)];
            labelmat_allsizes(newlabelmat>1) = newlabelmat(newlabelmat>1);
            blobs_labels = [blobs_labels; unique(newlabelmat(newlabelmat>1))];
        end
        for i = 1:length(blobs_labels)
            blobsizes(i) = length(find(labelmat_allsizes==blobs_labels(i)));
        end
        maxblobsizefound = max(blobsizes);
    end
    close(barhandle)
else
    rmfield(pars,{'visualize_shrinking', 'raisethresh_increment'});
end 

% eliminate patches smaller than minAreaPatch
if pars.minAreaPatch_squm > 0
    count = 0; 
    count2 = 0;
    patchimage_labeled = zeros(size(labelmat_allsizes));
    patchimage_sub_size_thresh = zeros(size(labelmat_allsizes));
    bnds_large = {};
    bnds_small = {};
    bnds_large_image = false(size(im));
    bnds_small_image = false(size(im));
    if ~exist('blobs_labels','var')
        blobs_labels = unique(labelmat_allsizes);
        blobs_labels(blobs_labels==0) = [];
    end
    for i = 1:length(bnds_allsizes)
        thisbloblabel = blobs_labels(i);
        blobnpix = length(find(labelmat_allsizes==thisbloblabel));
        if blobnpix >= minAreaPatch_pix
            count = count+1;
            bnds_large{count} = bnds_allsizes{i};
            labelval = labelmat_allsizes(bnds_allsizes{i}(1,1),bnds_allsizes{i}(1,2));
            inds = labelmat_allsizes == labelval;
            patchimage_labeled(inds) = labelval;
            bnds_large_image(sub2ind(size(im),bnds_allsizes{i}(:,1),bnds_allsizes{i}(:,2))) = true;
        else
            count2 = count+1;
            bnds_small{count2} = bnds_allsizes{i};
            inds = sub2ind(size(labelmat_allsizes),bnds_allsizes{i});
            labelval = labelmat_allsizes(bnds_allsizes{i}(1,1),bnds_allsizes{i}(1,2));
            inds = labelmat_allsizes == labelval;
            patchimage_sub_size_thresh(inds) = labelval;
            bnds_small_image(sub2ind(size(im),bnds_allsizes{i}(:,1),bnds_allsizes{i}(:,2))) = true;
        end
    end
    patchimage = logical(patchimage_labeled);
%     patchdata.meanPatchAreaSqum = patchdata.patches_area_sqmm * 1000^2 / length(bnds_large);
    patchdata.patchimage_sub_sizethresh = patchimage_sub_size_thresh;
    patchdata.bnds_sub_sizethresh = bnds_small_image;
    patchdata.bnds_patches = bnds_large_image;
    patchdata.patch_bnds = bnds_large;
    patchdata.sub_sizethresh_bnds = bnds_small;
elseif pars.minAreaPatch_squm <= 0
    patchdata.bnds_patches = edge(patchimage);
end
    
npixroi = length(find(roi));
patchdata.roi_area_sqmm = npixroi * [umPerPix(zoom,scope)]^2 / 1000^2; 
npixpatches = length(find(patchimage));
patchdata.patches_area_sqmm = npixpatches  * [umPerPix(zoom,scope)]^2 / 1000^2;
patchdata.fracAreaPatches = npixpatches / npixroi;
patchdata.patchimage = patchimage;
patchdata.nonpatchimage = roi & ~patchdata.patchimage; 

%% output results
patchdata.im = im;
patchdata.imfile = imfile;
patchdata.imroiblurred = imroiblurred; 
patchdata.imblurred = imblurred;
patchdata.minInRoi = minInRoi;
patchdata.maxInRoi = maxInRoi;
patchdata.pixPerUm = pixPerUm(zoom,scope);
patchdata.zoom = zoom;
patchdata.scope = scope;
patchdata.pars = pars;
if pars.save_patchdata
    save([save_prepend, 'patchdata'], 'patchdata')
    save_quant_borders_image(patchdata, pars.quantborders_ops)
end

%% plotting
if pars.show_plots
%     close all
%     figure('units','normalized','outerposition',[0 0 1 1])
%     imagesc(2*double(logical(labelmat_large)) + 1*double(logical(labelmat_small)))
%     colorbar
%     figure('units','normalized','outerposition',[0 0 1 1])
%     imagesc(imroiblurred)
%     colorbar
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    if strcmp(pars.threshMode,'intensityQuantiles')
       imroi = im;
       imroi(logical(quant_levels_img==0)) = 0; 
    end
    mm = max(double(imroi(:)));
    imagesc(double(imroi) + 2*mm*double(patchimage))
    colorbar
end
 
 
 