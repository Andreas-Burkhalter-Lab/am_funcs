%%% compare cell density from cellsFile in patches to interpatches
%%% from input patchData
%     patchResults = analyzePatchiness_cellDensity(patchAreaImageFile, cellsFile, zoom, intptchThicknessUM, patchImageFile)
% in:   patchAreaImageFile = image with white areas indicating patches, other areas black  
%       cellsFile = BW image from amygdala to POR, with m2... or can be the image matrix read from an image file  
%       zoom = zoom level of all input images
%       intptchThicknessUM = specified thickness of interpatches around drawn patches in microns
%       patchImageFile = m2 expression image (just for visual reference, not any processing) 
% % out:  patchResults = data structure containing areas of patches and
% %           interpatches, optical density results for patches and interpatches
% make sure input image is the same size and zoom as the data from which
%   patchData was taken
%%% To avoid nonzero pixel artifacts along the edge of area-defining BW images, set
%%% anti-aliasing in Illustrator to 'none'.
% last updated 18-05-09 on thermaltake

function patchResults = analyzePatchiness_cellDensity(patchAreaImageFile, cellsFile, zoom, intptchThicknessUM, patchImageFile)

if ~isempty(patchImageFile)
    patchImage = imread(patchImageFile);
    patchImage = patchImage(:,:,1); % in case image is rgb, take only first channel
else 
    patchImage = [];
end
[patchAreaImage patchAreaImageFile] = loadbw(cellsFile);
[cellsImage cellsFile] = loadbw(cellsFile);
xyPixelsPerUM = pixPerUm(zoom); 
sqMMPerPixel = xyPixelsPerUM^-2 / 10e6; % area of square pixel in sq mm
inptchThicknessPix = xyPixelsPerUM * intptchThicknessUM;
% get interpatch areas by outlining patch areas
[interpatchPixels_yx,patchPixels_yx,patchOutlines_yx,labelmat] = drawOutline(patchAreaImage,inptchThicknessPix,'noOverlap'); %#ok<ASGLU> % contains slow step

% check whether any area images have white pixels at borders, potentially
% indicating image processing artifacts
if any(find(patchAreaImage(:,1))) || any(find(patchAreaImage(:,end))) || any(find(patchAreaImage(1,:))) || any(find(patchAreaImage(end,:)))
    go_on = input(['Warning: patch area image ' patchAreaImageFile ' has noDnzero pixels along a border. Enter ''y'' to continue.'],'s');
    if ~strcmp(go_on,'y')
        error('quitting analyzePatchiness_cellDensity')
    end
end

% find cell centers
[bounds, labelmat, ncells_total, adjc] = bwboundaries(cellsimage);
cellcentersYX = NaN(ncells_total,2);
for j = 1:ncells_total
    [y x] = find(labelmat == j);
    cellcentersYX(j,1) = round(mean(y));
    cellcentersYX(j,2) = round(mean(x));
end
roiInds = find(roi);
cellCentersInds = sub2ind(size(roi),cellcentersYX(:,1),cellcentersYX(:,2));

npatches = length(interpatchPixels_yx); 
patchTable = table;
patchTable.patch_cellsPerSqMM = NaN(npatches,1); 
patchTable.interpatch_cellsPerSqMM = NaN(npatches,1);
patchTable.ratio_patchToIntptch = NaN(npatches,1);
patchTable.dif_patchToIntptch = NaN(npatches,1);
patchTable.nPatchCells = NaN(npatches,1); 
patchTable.nInterpatchCells = NaN(npatches,1);
patchTable.patchPixels_yx = patchPixels_yx; 
patchTable.interpatchPixels_yx = interpatchPixels_yx; 
patchTable.patchOutlines_yx = patchOutlines_yx; 

%% count cells

for patchind = 1:npatches
    %patches
    patchSubsYX = patchTable.patchPixels_yx{patchind};
    patchIndices = sub2ind(size(cellsImage),patchSubsYX(:,1),patchSubsYX(:,2));
    patchTable.nPatchCells(patchind) = length(intersect(patchIndices,cellCentersInds)); % number of cells in this patch
    nPatchPixels = length(patchIndices);
    patchAreaSqMM = nPatchPixels * sqMMPerPixel;
    patchTable.nPatchCells(patchind) = length(intersect(patchIndices,cellInds)); % number of cells in this patch
    patchTable.patch_cellsPerSqMM(patchind) = patchTable.nPatchCells(patchind) / patchAreaSqMM;
    %interpatches
    intpatchSubsYX = patchTable.interpatchPixels_yx{patchind};
    intpatchIndices = sub2ind(size(cellsImage),intpatchSubsYX(:,1),intpatchSubsYX(:,2));
    patchTable.nInterpatchCells(patchind) = length(intersect(intpatchIndices,cellCentersInds)); % number of cells in this interpatch
    nInterpatchPixels = length(intpatchIndices);
    interpatchAreaSqMM = nInterpatchPixels * sqMMPerPixel;
    patchTable.nInterpatchCells(patchind) = length(intersect(intpatchIndices,cellInds)); % number of cells in this patch
    patchTable.interpatch_cellsPerSqMM(patchind) = patchTable.nInterpatchCells(patchind) / interpatchAreaSqMM;
    %compare
    patchTable.ratio_patchToIntptch(patchind) = patchTable.patch_cellsPerSqMM(patchind) / patchTable.interpatch_cellsPerSqMM(patchind); % density ratio
    patchTable.dif_patchToIntptch(patchind) = patchTable.patch_cellsPerSqMM(patchind) - patchTable.interpatch_cellsPerSqMM(patchind); % density difference
end
    
patchResults.patchTable = patchTable;
patchResults.patchImage = patchImage;
patchResults.cellSubsYX = cellSubsYX;
patchResults.intptchThicknessUM = intptchThicknessUM;
patchResults.zoom = zoom;
patchResults.xyPixelsPerUM = xyPixelsPerUM;
patchResults.npatches = npatches; 
