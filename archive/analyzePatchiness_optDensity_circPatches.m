%% this version assumes that patches will be circular and interpatches will be circular rings around them
%%% compare intensity of input image imageFile in patches to interpatches
%%% from input patchData
%%%%% imageFile may be e.g. aav.egfp projections from amygdala to POR, with m2
%%%%% expression image in patchData
% make sure input image is the same size and zoom as the data from which
% patchData was taken
% last updated 1/30/17

function patchResults = analyzePatchiness_optDensity_circPatches(imageFile, patchData)

cmprsImage = imread(imageFile);
patchTable = patchData.patchTable;
npatches = height(patchTable);
patchTable.patchIntens = NaN(npatches,1);
patchTable.interpatchIntens = NaN(npatches,1);
xgrid = repmat(1:size(cmprsImage,2),size(cmprsImage,1),1);
ygrid = repmat([1:size(cmprsImage,1)]',1,size(cmprsImage,2));

for patchind = 1:npatches
    thisX = patchTable.centerx(patchind);
    thisY = patchTable.centery(patchind);
	thesePatchPixels = (xgrid-thisX).^2 + (ygrid-thisY).^2 < patchData.patchRadiusPix.^2;
	interpatchEnclosed = (xgrid-thisX).^2 + (ygrid-thisY).^2 < patchData.interpatchRadiusPix.^2;
    theseInterpatchPixels = interpatchEnclosed & ~thesePatchPixels;
    patchTable.patchIntens(patchind) = mean(cmprsImage(thesePatchPixels));
    patchTable.interpatchIntens(patchind) = mean(cmprsImage(theseInterpatchPixels));   
end

[patchData.h patchData.p] = ttest2(patchTable.patchIntens(:), patchTable.interpatchIntens(:));
patchData.patchTable = patchTable;
patchResults = patchData; 