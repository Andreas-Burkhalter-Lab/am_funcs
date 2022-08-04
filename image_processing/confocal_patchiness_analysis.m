%%% confocal_patchiness_analysis: get optical density of input to patches
%%% vs. interpatches in patch regions projected through a confocal stack
% in:   1. stackFile = .oib confocal stack file containing projection intensities
%       2. stackIndex = index of the stack to be analyzed within stackFile (eg. EGFP vs m2tdt stack) 
%       3. patchAreaImageFile = image with white areas indicating patches, other areas black  
%       4. areaBoundaryFile = BW image where white covers the entire area of interest (eg POR), the background is black
%       5. baselineAreaFile = BW image where white covers a region expected to have basline/zero fluorescence (eg. V1 for amyg feedback)   
%       6. zoom = zoom level of all input images on AB Lab epifluorescence scope
%       7. intptchThicknessUM = specified thickness of interpatches around drawn patches in microns
%%% updated 3/20/17 

function [stack_results, stack_info] = confocal_patchiness_analysis(stackFile, stackIndex, patchAreaImageFile, areaBoundaryFile, baselineAreaFile, zoom, intptchThicknessUM)

stackTable = getConfocalData(stackFile);
nslices = height(stackTable);
stack_info = struct;
stack_info.nslices = nslices;
stack_info.stackFile = stackFile; 
stack_info.stackIndex = stackIndex; 
stack_info.patchAreaImageFile = patchAreaImageFile; 
stack_info.areaBoundaryFile = areaBoundaryFile; 
stack_info.zoom = zoom; 
stack_info.intptchThicknessUM = intptchThicknessUM; 
patchres = cell(nslices,1);
pval = NaN(nslices,1);
mrat = NaN(nslices,1);
mdif = NaN(nslices,1);
stack_results = table(patchres,pval,mrat,mdif,'VariableNames',{'patch_results','pval','mratio_patchToIntptch','mdif_patchToIntptch'}); 

for sliceind = 1:nslices 
    sliceind
    sliceimg = stackTable{sliceind,stackIndex}{:};
    patchResults = analyzePatchiness_optDensity([], patchAreaImageFile, sliceimg, areaBoundaryFile, baselineAreaFile, zoom, intptchThicknessUM); 
    patchResults = rmfield(patchResults,{'patchImage','patchAreaImage','cmprsImage'});
    stack_results.patch_results{sliceind} = patchResults; 
    [~,stack_results.pval(sliceind)] = ttest(patchResults.patchTable.patchIntens, patchResults.patchTable.interpatchIntens); 
    stack_results.mratio_patchToIntptch(sliceind) = geomean(patchResults.patchTable.ratio_patchToIntptch);
    stack_results.mdif_patchToIntptch(sliceind) = mean(patchResults.patchTable.dif_patchToIntptch);
end

stack_info.scaleToPatchesIntptchs = patchResults.scaleToPatchesIntptchs;  
stack_info.npatches = patchResults.npatches;
