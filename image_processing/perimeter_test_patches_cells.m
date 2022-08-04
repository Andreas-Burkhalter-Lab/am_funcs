%%% compare cell density in patches to interpatches
%%% from input patchData
%    patchResults = analyzePatchiness_optDensity(m2ImageFile, patchAreaImageFile compareImageFile, zoom, pars.intptchThicknessUM)
% in:   patchAreaImageFile = image with white areas indicating patches, other areas black  
%       cellsFile = bw image containing one white patch for every cell to be counted 
%       areaBoundaryFile = BW image where white covers the entire area of interest (eg POR or V1), the background is black
%       zoom = zoom level of all input images
%       scope = scope used in acquiring images
%       pars = analysis parameters, including: normMethod, intptchThicknessUM, min_area_patch_sqmm
% out:  patchResults = data structure containing areas of patches and
%           interpatches, cell density results for patches and interpatches
% make sure input image is the same size and zoom as the data from which
% patchData was taken
%%% To avoid nonzero pixel artifacts along the edge of area-defining BW images, set
%%% anti-aliasing in Illustrator to 'none'.
% last updated 5/31/18

function patchResults = perimeter_test_patches_cells(patchimage, cellsFile, roi, zoom, scope, pars)
pars = vardefault('pars',struct);

pars.intptchThicknessUM = field_default(pars,'intptchThicknessUM',30); 
pars.min_area_patch_sqmm = field_default(pars,'min_area_patch_sqmm',0);

[patchimage, patchimageFile] = loadbw(patchimage); 
[roi, areaBoundaryImageFile] = loadbw(roi);
if ischar(cellsFile)
    cellsImage = imread(cellsFile);
    cellsImage = cellsImage(:,:,1); % in case image is rgb, take only first channel
else
    cellsImage = cellsFile;
end
if ~exist('scope','var')
    scope='epimicro';
end


% find cell centers
[bounds, labelmat, ncells_total, adjc] = bwboundaries(cellsImage);
cellCentersYX_raw = NaN(ncells_total,2);
for j = 1:ncells_total
    [y x] = find(labelmat == j);
    cellCentersYX_raw(j,1) = round(mean(y));
    cellCentersYX_raw(j,2) = round(mean(x));
end
cellCentersInds_raw = sub2ind(size(cellsImage),cellCentersYX_raw(:,1),cellCentersYX_raw(:,2));
%%% remove cells outside of the roi
roiinds = find(roi);
cellCentersInds = intersect(cellCentersInds_raw, roiinds);
[cellCentersYX(:,1), cellCentersYX(:,2)]  = ind2sub(size(cellsImage),cellCentersInds);

xyPixelsPerUM = pixPerUm(zoom,scope); 
sqMMPerPixel = xyPixelsPerUM^-2 / 10e6; % area of square pixel in sq mm
inptchThicknessPix = xyPixelsPerUM * pars.intptchThicknessUM;
minAreaPatch_pix = pars.min_area_patch_sqmm * (pixPerUm(zoom,scope))^2;
% get interpatch areas by outlining patch areas
[interpatchPixels_yx,patchPixels_yx,patchOutlines_yx,labelmat] = drawOutline(patchimage,inptchThicknessPix,'noOverlap'); %#ok<ASGLU> 

npatches = length(interpatchPixels_yx); 
patchTable = table;
patchTable.nPatchCells = NaN(npatches,1); 
patchTable.nInterpatchCells = NaN(npatches,1);
patchTable.ratio_patchToIntptch = NaN(npatches,1);
patchTable.patchPixels_yx = patchPixels_yx; 
patchTable.interpatchPixels_yx = interpatchPixels_yx; 
patchTable.patchOutlines_yx = patchOutlines_yx; 

for patchind = 1:npatches
        %patches
    patchSubsYX = patchTable.patchPixels_yx{patchind};
    patchIndices = sub2ind(size(cellsImage),patchSubsYX(:,1),patchSubsYX(:,2));
    patchTable.nPatchCells(patchind) = length(intersect(patchIndices,cellCentersInds)); % number of cells in this patch
    nPatchPixels = length(patchIndices);
    patchAreaSqMM = nPatchPixels * sqMMPerPixel;
    patchTable.patch_cellsPerSqMM(patchind) = patchTable.nPatchCells(patchind) / patchAreaSqMM;
    %interpatches
    intpatchSubsYX = patchTable.interpatchPixels_yx{patchind};
    intpatchIndices = sub2ind(size(cellsImage),intpatchSubsYX(:,1),intpatchSubsYX(:,2));
    patchTable.nInterpatchCells(patchind) = length(intersect(intpatchIndices,cellCentersInds)); % number of cells in this interpatch
    nInterpatchPixels = length(intpatchIndices);
    interpatchAreaSqMM = nInterpatchPixels * sqMMPerPixel;
    patchTable.interpatch_cellsPerSqMM(patchind) = patchTable.nInterpatchCells(patchind) / interpatchAreaSqMM;
    %compare
    patchTable.ratio_patchToIntptch(patchind) = patchTable.patch_cellsPerSqMM(patchind) / patchTable.interpatch_cellsPerSqMM(patchind); % density ratio
end
% don't include patches below the size threshold
patchTable.patchSizePix = cellfun(@(x)size(x,1), patchTable.patchPixels_yx);
deleterow = patchTable.patchSizePix <minAreaPatch_pix;
patchTable(deleterow,:) = [];
npatches = height(patchTable); % update after deleting small patches

     
[~, patchResults.pval_perimeter_test] = ttest(patchTable.patch_cellsPerSqMM, patchTable.interpatch_cellsPerSqMM); 
patchResults.patch_cellsPerSqMM_mean = mean(patchTable.patch_cellsPerSqMM);
patchResults.interpatch_cellsPerSqMM_mean = mean(patchTable.interpatch_cellsPerSqMM);
patchResults.patchInterpatchRatio = patchResults.patch_cellsPerSqMM_mean / patchResults.interpatch_cellsPerSqMM_mean;
patchResults.patchTable = patchTable;
patchResults.patchAreaImage = patchimage;
patchResults.cellsImage = cellsImage;
patchResults.intptchThicknessUM = pars.intptchThicknessUM;
patchResults.inptchThicknessPix = inptchThicknessPix;
patchResults.minAreaPatch_pix = minAreaPatch_pix;
patchResults.zoom = zoom;
patchResults.xyPixelsPerUM = xyPixelsPerUM;
patchResults.npatches = npatches; 
