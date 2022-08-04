%%% manually select patch centers from an input image 'patchImageFile' containing patches; this
%%% function will find the nearest local peak of intensity; will also
%%% output areas considered patches as 'patchData.patchPixels' (pixels
%%% within a radius of each patch center)
% 'patchImageData' = either an image file readable with imread.m  or (for confocal) 
%       a n-by-1 cell array of imported images of the same size (as from
%       getConfocalData.m); if there are multiple images, a composite of the stack will be usd
% 'xyPixelsPerUm' = pix/micron scale factor
%    -for confocal stacks, find by opening the
%    file with FIJI bioformats, divide '[Axis 0 Parameters Common] MaxSize'
%    by '[Axis 0 Parameters Common] EndPosition
%    -for sectioned tissue imaged on AB lab scope, use zoomConvert.m
% 'patchRadiusUM' selects the radius of patches in microns
% 'interpatchRadiusUM' selects the radius of interpatches in microns 
%       (from enclosed patch center to outer rim of enclosing interpatch) 
% 'confocalSlice' (only if patchImageFile is confocal .oib stack) picks
%       which slice or slices, counting from the top, will be selected from 
%       which to click on patches; if multiple slices are selected, a composite
%       image will be used
%%%%% last updated 1/30/17

function patchData = getPatches(patchImageData,xyPixelsPerUM,patchRadiusUM,interpatchRadiusUM)

patchRadiusPix = patchRadiusUM * xyPixelsPerUM;
interpatchRadiusPix = interpatchRadiusUM * xyPixelsPerUM;
outlineThresholdPatch = patchRadiusPix * 1.05; % for plotting the patch and interpatch outlines
outlineThresholdInterpatch = interpatchRadiusPix * 1.05; % for plotting the patch and interpatch outlines

%%% select patch centers manually
if iscell(patchImageData)
    sliceStack = inf(size(patchImageData{1},1),size(patchImageData{1},2),length(patchImageData));
    for sliceind = 1:length(patchImageData) % convert to 3d matrix
        sliceStack(:,:,sliceind) = patchImageData{sliceind};
    end
    patchImage = mean(sliceStack,3); % make composite
elseif ischar(patchImageData)
    patchImage = imread(patchImageData);
end
close all
figr = figure;
imagesc(patchImage)
hold on
fprintf('Left-click on patch centers, then right-click when all patch centers are selected.\n')

clickedPoint = [NaN NaN]; % init
clickedPatchCenters = []; % init
buttonPressed = 1; % init
while buttonPressed == 1; % press Enter to make clickedPoint empty
    [clickedPoint(1) clickedPoint(2) buttonPressed] = ginput(1);
    if buttonPressed == 1
        clickedPatchCenters = [clickedPatchCenters; clickedPoint]; % add point to list
        scatter(clickedPatchCenters(:,1),clickedPatchCenters(:,2),'r')
    end
end

%%% get patch areas
% use manual patch centers for now
npatches = size(clickedPatchCenters,1);
patchPixels = false(size(patchImage));
patchOutlines = false(size(patchImage));
interpatchPixels = false(size(patchImage));
interpatchOutlines = false(size(patchImage));
xgrid = repmat(1:size(patchImage,2),size(patchImage,1),1);
ygrid = repmat([1:size(patchImage,1)]',1,size(patchImage,2));
for patchInd = 1:npatches
    thisX = clickedPatchCenters(patchInd,1);
    thisY = clickedPatchCenters(patchInd,2);
    thesePatchPixels = (xgrid-thisX).^2 + (ygrid-thisY).^2 < patchRadiusPix.^2;
    thisPatchOutline = abs((xgrid-thisX).^2 + (ygrid-thisY).^2 - patchRadiusPix.^2) < outlineThresholdPatch;
    theseInterpatchPixels = (xgrid-thisX).^2 + (ygrid-thisY).^2 < interpatchRadiusPix.^2 ...
                            & ~thesePatchPixels; 
    thisInterpatchOutline = abs((xgrid-thisX).^2 + (ygrid-thisY).^2 - patchRadiusPix.^2) < outlineThresholdInterpatch;
    patchPixels = patchPixels | thesePatchPixels;
    patchOutlines = patchOutlines | thisPatchOutline;
    interpatchPixels = interpatchPixels | theseInterpatchPixels;
    interpatchOutlines = interpatchOutlines | thisInterpatchOutline;
end
hold off
close(figr)

patchTable = table(clickedPatchCenters(:,1),clickedPatchCenters(:,2),'VariableNames',{'centerx','centery'});
patchData.patchTable = patchTable;
patchData.patchPixels = patchPixels;
patchData.interpatchPixels = interpatchPixels;
patchData.xyPixelsPerUM = xyPixelsPerUM; 
patchData.patchRadiusUM = patchRadiusUM;
patchData.interpatchRadiusUM = interpatchRadiusUM;
patchData.patchRadiusPix = patchRadiusPix;
patchData.interpatchRadiusPix = interpatchRadiusPix;
patchData.patchImageData = patchImageData; 
patchData.patchImage = patchImage;



    