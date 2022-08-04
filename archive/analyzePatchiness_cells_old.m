%  patchResults = analyzePatchiness_cells(cellStack, patchData)
%%% compare cell count density of cellStack in patches to interpatches
%%% from input patchData - all cells underneath a patch/interpatch will be
%%% counted as being in the patch/interpatch
%%%% in 1: cellStack = x-by-1 cell array with each row containing cell
%%%%    coordinates from each slice (cellsTable.goodCells column from removeDuplicateClickedCells.m)  
%%%% in 2: patchData output from getpatches
% last updated 12/6/16

function patchResults = analyzePatchiness_cells_old(cellStack, patchData)

patchTable = patchData.patchTable;
patchImage = patchData.patchImage;
npatches = height(patchTable);
patchTable.patchCellsPerUMSq = NaN(npatches,1); 
patchTable.interpatchCellsPerUMSq = NaN(npatches,1);
patchAreaUMSq = pi * patchData.patchRadiusUM^2;
interpatchEnclosedAreaUMSq = pi * patchData.interpatchRadiusUM^2;
interpatchAreaUMSq = interpatchEnclosedAreaUMSq - patchAreaUMSq;
xgrid = repmat(1:size(patchImage,2),size(patchImage,1),1);
ygrid = repmat([1:size(patchImage,1)]',1,size(patchImage,2));
allCellsUM = cell2mat(cellStack); %unpack - xyz coordinates in microns
allCellsPix = round(patchData.xyPixelsPerUM * allCellsUM); % coordinates in pixels
ncells = size(allCellsPix,1);

% count cells in patches and interpatches
for patchind = 1:npatches
    thisX = patchTable.centerx(patchind);
    thisY = patchTable.centery(patchind);
	thesePatchPixels = (xgrid-thisX).^2 + (ygrid-thisY).^2 < patchData.patchRadiusPix.^2;
    [yPatch xPatch] = find(thesePatchPixels); % patch pixels as coordinates
    thesePatchPixels_xy = [xPatch yPatch];
	interpatchEnclosed = (xgrid-thisX).^2 + (ygrid-thisY).^2 < patchData.interpatchRadiusPix.^2;
    theseInterpatchPixels = interpatchEnclosed & ~thesePatchPixels;
    [yInterpatch xInterpatch] = find(theseInterpatchPixels); % patch pixels as coordinates
    theseInterpatchPixels_xy = [xInterpatch yInterpatch];
    cellIsWithinPatch = false(ncells,1); % turn to true if the cell in this index is within the patch
    cellIsWithinInterpatch = false(ncells,1); % turn to true if the cell in this index is within the interpatch
    for cellInd = 1:ncells
        xMatchThisPatch = allCellsUM(cellInd,1) == thesePatchPixels_xy(:,1);
        yMatchThisPatch = allCellsUM(cellInd,2) == thesePatchPixels_xy(:,2);
        cellIsWithinPatch(cellInd) = any(xMatchThisPatch & yMatchThisPatch); % check if cell coord matches any patch coord
        xMatchThisInterpatch = allCellsUM(cellInd,1) == theseInterpatchPixels_xy(:,1);        
        yMatchThisInterpatch = allCellsUM(cellInd,2) == theseInterpatchPixels_xy(:,2); 
        cellIsWithinPatch(cellInd) = any(xMatchThisInterpatch & yMatchThisInterpatch);% check if cell coord matches any interpatch coord
    end
    nCellsInPatch = length(find(cellIsWithinPatch));
    nCellsInInterpatch = length(find(cellIsWithinInterpatch));  
    patchTable.patchCellsPerUMSq(patchind) = nCellsInPatch / patchAreaUMSq; % cell density in um squared
    patchTable.interpatchCellsPerUMSq(patchind) = nCellsInInterpatch / interpatchAreaUMSq;  % cell density in um squared
    clear cellIsWithinPatch cellIsWithinInterpatch nCellsInPatch nCellsInInterpatch
end

[patchData.h patchData.p] = ttest2(patchTable.patchCellsPerUMSq(:), patchTable.interpatchCellsPerUMSq(:));
patchData.patchTable = patchTable;
patchResults = patchData;