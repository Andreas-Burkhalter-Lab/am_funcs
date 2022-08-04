%%% analyze patchiness of optical density from multiple image sets and
%%% compile the results
% updated 3/10/17
%%% inputs: 
%%%     filelist = table with columns labeled 'subject,m2_image,feedback_image,patch_borders_image,roi_area_image,baseline_area_image,zoom'  
%%%     intptchThicknessUM = thickness of interpatches to draw around patches in micron   

function [batch_results batch_info] = batch_patchiness_analysis(filelist,intptchThicknessUM)

nsubs = height(filelist);
batch_results = table; % concat patches from subjects
batch_info = struct;

for subind = 1:nsubs
    patchResults = analyzePatchiness_optDensity(filelist.m2_image{subind},filelist.patch_borders_image{subind}, ...
         filelist.feedback_image{subind}, filelist.roi_area_image{subind}, filelist.baseline_area_image{subind}, filelist.zoom(subind), intptchThicknessUM);
    patchTable = patchResults.patchTable; 
    patchTable.subject = repmat(filelist.subject(subind),patchResults.npatches,1);
    batch_results = [batch_results; patchTable]; %#ok<AGROW> % concat subjects
    batch_info(subind) = patchResults;
end