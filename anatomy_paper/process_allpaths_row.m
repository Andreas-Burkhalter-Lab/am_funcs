%%%% run patch finding and patch vs interpatch comparison for entry in
%%%% allpaths table
% updated 2020/3/10

row_to_process = 6;
sort_rows = 0; % sort rows alphabetically by pathway

for thisvar = {'m2_orig' 'proj_orig' 'm2filt' 'projfilt' 'roi' 'baseline'}
    thisvar = thisvar{:};
    varval = allpaths{row_to_process,thisvar}{:};
    if ~isempty(varval) && isempty(fileparts(varval)) % add path if needed
        allpaths{row_to_process,thisvar}{:} = [pwd filesep varval];
    end
end

%%%% following line may cause error when commented in... use findpatches.m instead
% allpaths.patches{row_to_process} = analyzePatches(allpaths.m2filt{row_to_process},allpaths.roi{row_to_process},allpaths.zoom{row_to_process},allpaths.scope{row_to_process});
allpaths.patches{row_to_process}= findpatches(allpaths.m2filt{row_to_process},allpaths.roi{row_to_process},allpaths.zoom(row_to_process),allpaths.scope{row_to_process}); % updated function

if ~isempty( allpaths.proj_orig{row_to_process} )
    if ~isfield(allpaths.patches{row_to_process},'labelmat_large') % add labelmat if it wasn't added during findpathes
        allpaths.patches{row_to_process}.labelmat_large = bwlabel(allpaths.patches{row_to_process}.patchimage,8);
    end
    allpaths.permtest_projorig{row_to_process,1} = permutation_test_patches(logical(allpaths.patches{row_to_process}.labelmat_large),allpaths.proj_orig{row_to_process},...
        allpaths.roi{row_to_process},allpaths.baseline{row_to_process},allpaths.zoom(row_to_process),allpaths.scope{row_to_process});
    allpaths.patchInterpatchRatio_projorig(row_to_process,1) = allpaths.permtest_projorig{row_to_process}.patchInterpatchRatio;
end
allpaths.permtest_projfilt{row_to_process,1} = permutation_test_patches(logical(allpaths.patches{row_to_process}.labelmat_large),allpaths.projfilt{row_to_process},...
    allpaths.roi{row_to_process},allpaths.baseline{row_to_process},allpaths.zoom(row_to_process),allpaths.scope{row_to_process});
allpaths.patchInterpatchRatio_projfilt(row_to_process,1) = allpaths.permtest_projfilt{row_to_process}.patchInterpatchRatio;
patchoutlines_save_name = strrep(allpaths.roi{row_to_process},'area.png','patchoutlines'); 
% savetif(allpaths.patches{row_to_process}.bnds_large_image,patchoutlines_save_name)

if sort_rows
    [~,rr] = sort(char(allpaths.path)); 
    allpaths = allpaths(rr(:,1),:);
end
    
% % % % save patch borders image
% save_quant_borders_image(allpaths.patches{row_to_process})

clear thisvar varval patchoutlines_save_name row_to_process