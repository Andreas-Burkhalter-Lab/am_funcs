 % compare cell density of in vs outside of patches and
 % perform permuation test
 %%% results = permutation_test_patches(patchFile, projFile, roiFile, baselineFile, zoom, scope, pars)
 % updated 18/6/14 
 function results = permutation_test_patches_cells(patchimage, cellsImage, roi, zoom, scope, pars)
 pars = vardefault('pars',struct);
 
 pars.npermutations = field_default(pars,'npermutations',2000); % number of permutations
    
[cellsImage, cellsImageFile] = loadbw(cellsImage);
results.cellsImage = cellsImage;
[patchimage, patchImageFile] = loadbw(patchimage);
results.patchimg = patchimage;
if isfield(pars,'interpatchimage')
    fprintf('Using custom interpatch image...')
    results.custom_interpatches = true;
    [interpatchimage, interpatchimage_file] = loadbw(pars.interpatchimage);
    results.interpatchimage = interpatchimage;
else
    results.custom_interpatches = false;
end
    
if exist('roi','var') && ~isempty(roi)
    [roi, roiFile] = loadbw(roi);
elseif ~exist('roi','var') || isempty(roi) % analyze all pixels
    results.roiFile = [];
    roi = true(size(cellsImage));
end
results.roi = roi;
    
if any(any(repmat(size(patchimage),2,1) ~= ([size(cellsImage); size(roi)])))
    error('image size mismatch')
end

xyPixelsPerUM = pixPerUm(zoom,scope); 
sqMMPerPixel = xyPixelsPerUM^-2 / 10e6; % area of square pixel in sq mm


%% permutation test
[cellCentersInds_raw, cellCentersInds_roiRelative, cellCentersImage] = find_shape_centers(cellsImage,roi); % find cell centers
patch_dummyvals = patchimage(roi); % get indices of patch pixels WITHIN the roi, as a logical vector
if isfield(pars,'interpatchimage') % custom interpatches
    interpatch_dummyvals = interpatchimage(roi); % get indices of custom interpatch pixels WITHIN the roi, as a logical vector
else % use all nonpatch roi pixels as interpatchpixels
    interpatch_dummyvals = ~patch_dummyvals; % get indices of interpatch pixels WITHIN the roi, as a logical vector
end 

n_patchpix = length(find(patch_dummyvals)); 
n_interpatchpix = length(find(interpatch_dummyvals)); 
patch_area_sqmm = n_patchpix * sqMMPerPixel;
interpatch_area_sqmm = n_interpatchpix * sqMMPerPixel;
npatchcells = nnz(patch_dummyvals & cellCentersInds_roiRelative); % total number of cells within roi patches
ninterpatchcells = nnz(interpatch_dummyvals & cellCentersInds_roiRelative); % total number of cells within roi interpatches
patch_cellsPerSqMM = npatchcells / patch_area_sqmm;
interpatch_cellsPerSqMM = ninterpatchcells / interpatch_area_sqmm;
patchInterpatchRatio = patch_cellsPerSqMM / interpatch_cellsPerSqMM; 
patchInterpatchDif = patch_cellsPerSqMM-interpatch_cellsPerSqMM;

shuffledist_patchInterpatchRatio = NaN(pars.npermutations,1);
wbar = waitbar(0,'Performing permutation test on cells in patches vs. interpatches...');
for i = 1:pars.npermutations %%% run the permutation test
    shuffled_cells = cellCentersInds_roiRelative(randperm(length(cellCentersInds_roiRelative))); 
    shuffled_npatchcells = nnz(patch_dummyvals & shuffled_cells); % total number of cells within roi patches
    shuffled_ninterpatchcells = nnz(interpatch_dummyvals & shuffled_cells); % total number of cells within roi interpatches
    shuffled_patch_cellsPerSqMM = shuffled_npatchcells / patch_area_sqmm;
    shuffled_interpatch_cellsPerSqMM = shuffled_ninterpatchcells / interpatch_area_sqmm;
    shuffledist_patchInterpatchRatio(i) = shuffled_patch_cellsPerSqMM / shuffled_interpatch_cellsPerSqMM; % get patch/interpatch density ratio for this iteration
    if any(i == round( [.01:.01:1] * pars.npermutations))
        try waitbar(i/pars.npermutations,wbar); end
    end
end
try close(wbar); end

%% outputs 
results.patch_cellsPerSqMM = patch_cellsPerSqMM;
results.interpatch_cellsPerSqMM = interpatch_cellsPerSqMM;
results.shuffledist_patchInterpatchRatio = shuffledist_patchInterpatchRatio;
results.shufflemean_patchInterpatchRatio = mean(shuffledist_patchInterpatchRatio);
results.patchInterpatchRatio = patchInterpatchRatio; 
results.patchInterpatchDif = patchInterpatchDif;
%%% two-tailed significance test: for each ratio, get abs(ratio-1) to get its deviation from 1
%%%     -use these values to get probability of finding ratio at least as different from 1 as the ratio that was found
results.pval_permutation_test = length(find(abs(shuffledist_patchInterpatchRatio-1)>abs(patchInterpatchRatio-1))) / pars.npermutations;
results.nperm = pars.npermutations;


% rearrange fields
if exist('movevars','file')
    firstfields = {'pval_permutation_test','patchInterpatchRatio','patchInterpatchDif','shufflemean_patchInterpatchRatio','nperm'};
    results = table2struct(movevars(struct2table(results,'AsArray',true), firstfields,'Before',1));
end





