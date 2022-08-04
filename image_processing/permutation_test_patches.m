 % compare intensity of projimage pixels in vs outside of patches and
 % perform permuation test
 %%% results = permutation_test_patches(patchFile, projFile, roiFile, baselineFile, zoom, scope, pars)
 % updated 2020/12/3
 function results = permutation_test_patches(patchFile, projFile, roiFile, baselineFile, zoom, scope, pars)
 
 pars = vardefault('pars',struct);
 pars.npermutations = field_default(pars,'npermutations',1000); % number of permutations
 pars.resample_pix = field_default(pars,'resample_pix',true); %%% resample pixels in proj image; set pixel area to force below
    % at epi scope x3.2, 1 pixel area == 3.84 squm.... area of a 60um-radius patch == 3600 squm
    pars.resampled_pix_area_squm = field_default(pars,'resampled_pix_area_squm',25); % area of a pixel in square microns after resampling
 pars.loadbw_pars = field_default(pars,'loadbw_pars',struct); 
 
%  pars.normMethod = field_default(pars,'none');
% pars.normMethod = field_default(pars,'divideByRoiMean');
pars.normMethod = field_default(pars,'normMethod','subtractBaseline_divideByRoiMean');

if ischar(projFile)
    projimg = imread(projFile);
    results.projFile = projFile;
else
    projimg = projFile;
end
projimg = projimg(:,:,1); % in case image is rgb, take only first channel
results.comparison_image = projimg;
 
if ischar(patchFile)
    patchimage = imread(patchFile);
    results.patchFile = patchFile;
else
    patchimage = patchFile; clear patchFile
end
patchimage = patchimage(:,:,1); % in case image is rgb, take only first channel
if numel(unique(patchimage(patchimage~=0))) ~= 1
    error('More or less than 1 unique non-zero pixel value found in patch file.')
end
patchimage = logical(patchimage); 
results.patchimg = patchimage;

if isfield(pars,'interpatchimage')
    fprintf('Using custom interpatch image...')
    results.custom_interpatches = true;
    [interpatchimage, interpatchimage_file] = loadbw(pars.interpatchimage);
    results.interpatchimage = interpatchimage;
else
    results.custom_interpatches = false;
end
    
if exist('roiFile','var') && ~isempty(roiFile)
    if ischar(roiFile)
        roi = imread(roiFile);
        results.roiFile = roiFile;
    else
        roi = roiFile; clear roiFile
    end
    
    roi = roi(:,:,1); % in case image is rgb, take only first channel
    if numel(unique(roi(roi~=0))) ~= 1
        error('More or less than 1 unique non-zero pixel value found in ROI file.')
    end
    roi = roi>0; 
elseif ~exist('roiFile','var') || isempty(roiFile) % analyze all pixels
    results.roiFile = [];
    roi = true(size(projimg));
end
results.roi = roi;

if exist('baselineFile','var') && ~isempty(baselineFile) % if baseline file provided
    baselineimg = loadbw(baselineFile, pars.loadbw_pars); 
    if any(size(patchimage) ~= size(baselineimg))
        error('image size mismatch')
    end
elseif strcmp(pars.normMethod, 'subtractBaseline_divideByRoiMean') % if baseline area is required but not provided, get roi manually
    fig1 = figure;
    imagesc(projimg)
    fprintf('Draw outline of baseline area. Double click to complete.')
    baselineimg = false(size(projimg));
    while isempty(find(baselineimg))
        baselineimg = roipoly;
        if isempty(find(baselineimg))
            fprintf('\nBaseline area too small, please draw larger.')
        end
    end
    close(fig1)
end
    
if any(any(repmat(size(patchimage),2,1) ~= ([size(projimg); size(roi)])))
    error('image size mismatch')
end

%% resample
results.resample_pix = pars.resample_pix;
if pars.resample_pix
    if ~exist('scope','var') || isempty(scope)
        scope = 'epimicro';
    end
    original_pix_area_squm = umPerPix(zoom,scope)^2;
    resample_factor = sqrt(original_pix_area_squm / pars.resampled_pix_area_squm);
    patchimage = imresize(patchimage,resample_factor);
    roi = imresize(roi,resample_factor);
    projimg = imresize(projimg,resample_factor);
    results.patchimg_resampled = patchimage;
    results.roi_resampled = roi;
    results.projimg_resampled = projimg;
    results.resampled_pix_area_squm = pars.resampled_pix_area_squm; 
    if strcmp(pars.normMethod, 'subtractBaseline_divideByRoiMean')
        results.baselineimg = baselineimg; 
        baselineimg = imresize(baselineimg,resample_factor);
        results.baselineimg_resampled = baselineimg; 
    end
    if isfield(pars,'interpatchimage')
        interpatchimage = imresize(interpatchimage, resample_factor);
        results.interpatchimage_resampled = interpatchimage;
    end
elseif ~pars.resample_pix
    results.baselineimg = baselineimg;     
end

%% norm
projroivals = double(projimg(roi));
switch pars.normMethod
    case 'none' 
        projroivals = projroivals;
    case 'divideByRoiMean'
        roiMean = mean(projroivals); 
        results.roiMean = roiMean;
        projroivals = (1/roiMean) * projroivals;  
    case 'subtractBaseline_divideByRoiMean'
        if ~isempty(find(baselineimg)) % if baseline exists after imresize
            basevals = projimg(baselineimg); % use resized baseline and resized projimg
        elseif isempty(find(baselineimg)) % if baseline was too small disappeared during imresize
            basevals = results.comparison_image(results.baselineimg); % use baseline from before imresize
        end

        baseline = mean(basevals);
        projroivals = projroivals - baseline;
        roiMean = mean(projroivals); 
        projroivals = (1/roiMean) * projroivals;  
        results.baseline = baseline;
        results.roiMean = roiMean;
end

%% perform permutation test
patch_dummyvals = patchimage(roi); % get indices of patch pixels WITHIN the roi, as a logical vector
if isfield(pars,'interpatchimage') % custom interpatches
    interpatch_dummyvals = interpatchimage(roi); % get indices of custom interpatch pixels WITHIN the roi, as a logical vector
else % use all nonpatch roi pixels as interpatchpixels
    interpatch_dummyvals = ~patch_dummyvals;
end
patchmean = mean(projroivals(patch_dummyvals));
interpatchmean = mean(projroivals(interpatch_dummyvals));
patchInterpatchDif = patchmean-interpatchmean;
patchInterpatchRatio = patchmean / interpatchmean; 
n_patchpix = length(find(patch_dummyvals)); 
n_interpatchpix = length(find(interpatch_dummyvals)); 
npixtotal = length(patch_dummyvals); % total number of pix in the roi

shuffledist_patchInterpatchRatio = NaN(pars.npermutations,1);
wbar = waitbar(0,'Performing permutation test on intensity in patches...');
for i = 1:pars.npermutations
    perm = randperm(npixtotal);
    newpatchpixinds = perm(1:n_patchpix);
    newinterpatchpixinds = perm(n_patchpix+1 : n_patchpix+n_interpatchpix); % note: if using custom interpatches, n_patchpix+n_interpatchpix ~= total roi pix
    newpatchmean = mean(projroivals(newpatchpixinds));
    newinterpatchmean = mean(projroivals(newinterpatchpixinds));
    shuffledist_patchInterpatchRatio(i) = newpatchmean / newinterpatchmean;
    if any(i == round( [.01:.01:1] * pars.npermutations))
        waitbar(i/pars.npermutations,wbar)
    end
end
close(wbar)

%% clean up
results.patchmean = patchmean;
results.interpatchmean = interpatchmean;
results.shuffledist_patchInterpatchRatio = shuffledist_patchInterpatchRatio;
results.shufflemean_patchInterpatchRatio = mean(shuffledist_patchInterpatchRatio);
results.patchInterpatchRatio = patchInterpatchRatio; 
results.patchInterpatchDif = patchInterpatchDif;
results.pval_permutation_test = length(find(abs(shuffledist_patchInterpatchRatio-1)>abs(patchInterpatchRatio-1))) / pars.npermutations;


results.nperm = pars.npermutations;
results.normMethod = pars.normMethod;

% rearrange fields
if exist('movevars','file')
    firstfields = {'pval_permutation_test','patchInterpatchRatio','patchInterpatchDif','shufflemean_patchInterpatchRatio','nperm'};
    results = table2struct(movevars(struct2table(results,'AsArray',true), firstfields,'Before',1));
end





