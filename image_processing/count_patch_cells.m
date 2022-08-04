%%% COUNT_PATCH_CELLS: count cells in patches vs interpatches in a
%%% list of labeled cell images
% [ftable_out] = countcells(filetable_file,patchdata_file)
%%% 

function [ftable_out] = count_patch_cells(filetable_file,patchdata_file)

load(patchdata_file); % contains patchdata
load(filetable_file); % contains filetable
ftable_out = filetable;
patchimage = patchdata.patchimage;
interpatchimage = patchdata.interpatchimage;


ftable_out = countcells(ftable_out,patchimage);
vars_to_rename = {'ncells','cellsPerSqMM','cellsYX'};
for varn = vars_to_rename
    ftable_out.Properties.VariableNames{strcmp(varn,ftable_out.Properties.VariableNames)} = ['patch_' varn{1}];
end
ftable_out = countcells(ftable_out,interpatchimage);
for varn = vars_to_rename
    ftable_out.Properties.VariableNames{strcmp(varn,ftable_out.Properties.VariableNames)} = ['ipatch_' varn{1}];
end