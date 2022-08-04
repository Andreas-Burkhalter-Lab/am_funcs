%%% compare intensity of input image imageFile in patches to interpatches
%%% from input patchData
%    patchResults = analyzePatchiness_optDensity(m2ImageFile, patchAreaImageFile compareImageFile, zoom, pars.intptchThicknessUM)
% in:   patchAreaImageFile = image with white areas indicating patches, other areas black  
%       compareImageFile = aav projections from amygdala to POR... or can be the image matrix read from an image file  
%       areaBoundaryFile = BW image where white covers the entire area of interest (eg POR or V1), the background is black
%       baselineAreaFile = BW image where white covers a region expected to have baseline/zero fluorescence (eg. V1 for amyg feedback)   
%       zoom = zoom level of all input images
%       scope = scope used in acquiring images
%       pars = analysis parameters, including: normMethod, intptchThicknessUM, min_area_patch_sqmm
% out:  patchResults = data structure containing areas of patches and
%           interpatches, optical density results for patches and interpatches
% make sure input image is the same size and zoom as the data from which
% patchData was taken
%%% To avoid nonzero pixel artifacts along the edge of area-defining BW images, set
%%% anti-aliasing in Illustrator to 'none'.
% last updated 5/23/18

function patchResults = perimeter_test_patches(patchimage, compareImageFile, areaBoundaryImage, baselineAreaImage, zoom, scope, pars)
pars = vardefault('pars',struct);

% pars.normMethod = field_default(pars,'normMethod','subtractBaseline_divideByRoiMean'); % scale to average intensity within roi
% pars.normMethod = field_default(pars,'normMethod','subtractBaseline_divideByRoiMax'); % scale to max intensity within roi
% pars.normMethod = field_default(pars,'normMethod','subtractBaseline_divideByPatchesIntptchs'); % scale intensity to avg of patches and interpatches, not avg of ROI
pars.normMethod = field_default(pars,'normMethod','subtractBaseline'); % use the raw, unscaled intensity values

pars.intptchThicknessUM = field_default(pars,'intptchThicknessUM',30); 
pars.min_area_patch_sqmm = field_default(pars,'min_area_patch_sqmm',0);

[patchimage, patchimageFile] = loadbw(patchimage); 
[areaBoundaryImage, areaBoundaryImageFile] = loadbw(areaBoundaryImage);
[baselineAreaImage, baselineAreaFile] = loadbw(baselineAreaImage);
if ischar(compareImageFile)
    cmprsImage = imread(compareImageFile);
    cmprsImage = cmprsImage(:,:,1); % in case image is rgb, take only first channel
else
    cmprsImage = compareImageFile;
end
if ~exist('scope','var')
    scope='epimicro';
end

xyPixelsPerUM = pixPerUm(zoom,scope); 
inptchThicknessPix = xyPixelsPerUM * pars.intptchThicknessUM;
% get interpatch areas by outlining patch areas
[interpatchPixels_yx,patchPixels_yx,patchOutlines_yx,labelmat] = drawOutline(patchimage,inptchThicknessPix,'noOverlap'); %#ok<ASGLU> % contains slow step
% get mean, max intensity across the entire area of interest
roiPixLocs = find(areaBoundaryImage);
roiPixVals = cmprsImage(roiPixLocs); %#ok<*FNDSB>
roiMean = mean(roiPixVals);
roiMax = max(roiPixVals);
% get mean intensity in baseline area
basePixLocs = find(baselineAreaImage);
basePixVals = cmprsImage(basePixLocs);
baseline = mean(basePixVals); 
minAreaPatch_pix = pars.min_area_patch_sqmm * (pixPerUm(zoom,scope))^2;

npatches = length(interpatchPixels_yx); 
patchTable = table;
patchTable.patchIntens = NaN(npatches,1); 
patchTable.interpatchIntens = NaN(npatches,1);
patchTable.ratio_patchToIntptch = NaN(npatches,1);
patchTable.dif_patchToIntptch = NaN(npatches,1);
patchTable.patchPixels_yx = patchPixels_yx; 
patchTable.interpatchPixels_yx = interpatchPixels_yx; 
patchTable.patchOutlines_yx = patchOutlines_yx; 

for patchind = 1:npatches
    patchSubsYX = patchTable.patchPixels_yx{patchind};
    patchIndices = sub2ind(size(cmprsImage),patchSubsYX(:,1),patchSubsYX(:,2));
    patchIntens_values_raw = cmprsImage(patchIndices); % all intensity values in the comparison image in this interpatch
    patchIntens_values_baseSubtracted = patchIntens_values_raw - baseline; % normalize to baseline
    if strcmp(pars.normMethod,'subtractBaseline_divideByRoiMean')
        patchIntens_values = patchIntens_values_baseSubtracted ./ roiMean; % get as multiple of average intensity in the ROI
    else 
        patchIntens_values = patchIntens_values_baseSubtracted;
    end
    patchTable.patchIntens(patchind) = mean(patchIntens_values); % normalized and scaled density of intensity in this interpatch
    intpatchSubsYX = patchTable.interpatchPixels_yx{patchind};
    intpatchIndices = sub2ind(size(cmprsImage),intpatchSubsYX(:,1),intpatchSubsYX(:,2));
    interpatchIntens_values_raw = cmprsImage(intpatchIndices); % all intensity values in the comparison image in this interpatch
    interpatchIntens_values_baseSubtracted = interpatchIntens_values_raw - baseline; % normalize to baseline
    if strcmp(pars.normMethod,'subtractBaseline_divideByRoiMean')
        interpatchIntens_values = interpatchIntens_values_baseSubtracted ./ roiMean; % get as multiple of average intensity in the ROI
    else
        interpatchIntens_values = interpatchIntens_values_baseSubtracted; 
    end
    patchTable.interpatchIntens(patchind) = mean(interpatchIntens_values); % density of intensity in this interpatch
    patchTable.ratio_patchToIntptch(patchind) = patchTable.patchIntens(patchind) / patchTable.interpatchIntens(patchind);
    patchTable.dif_patchToIntptch(patchind) = patchTable.patchIntens(patchind) - patchTable.interpatchIntens(patchind);
end
patchTable.patchSizePix = cellfun(@(x)size(x,1), patchTable.patchPixels_yx);
deleterow = patchTable.patchSizePix <minAreaPatch_pix;
patchTable(deleterow,:) = [];

switch pars.normMethod
    case 'subtractBaseline_divideByRoiMax'
        patchTable.patchIntens = patchTable.patchIntens ./ double(roiMax); 
        patchTable.interpatchIntens = patchTable.interpatchIntens ./ double(roiMax);
        patchTable.dif_patchToIntptch = patchTable.dif_patchToIntptch ./ double(roiMax);
    case 'subtractBaseline_divideByPatchesIntptchs'
        patchIntptchMean = mean([patchTable.patchIntens; patchTable.interpatchIntens]);
        patchTable.patchIntens = patchTable.patchIntens ./ patchIntptchMean; 
        patchTable.interpatchIntens = patchTable.interpatchIntens ./ patchIntptchMean;
        patchTable.dif_patchToIntptch = patchTable.dif_patchToIntptch ./ patchIntptchMean;
end
     
[~, patchResults.pval_perimeter_test] = ttest(patchTable.patchIntens, patchTable.interpatchIntens); 
patchResults.patch_mean = mean(patchTable.patchIntens);
patchResults.interpatch_mean = mean(patchTable.interpatchIntens);
patchResults.patchInterpatchRatio = patchResults.patch_mean / patchResults.interpatch_mean;
patchResults.patchTable = patchTable;
patchResults.patchAreaImage = patchimage;
patchResults.cmprsImage = cmprsImage;
patchResults.intptchThicknessUM = pars.intptchThicknessUM;
patchResults.inptchThicknessPix = inptchThicknessPix;
patchResults.normMethod = pars.normMethod; 
patchResults.minAreaPatch_pix = minAreaPatch_pix;
patchResults.zoom = zoom;
patchResults.xyPixelsPerUM = xyPixelsPerUM;
patchResults.npatches = npatches; 
patchResults.baseline = vardefault('baseline',0);
