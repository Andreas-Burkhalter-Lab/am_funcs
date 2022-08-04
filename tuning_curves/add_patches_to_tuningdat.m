%%% add patchdata to already-processed tuning data that doesn't have patch locations yet
% input the excel analysis master file name and desired rows, as with batch_full_session_analysis.m
%
%  add_patches_to_tuningdat(filetable_in, excel_rows)
%
% updated 2020/8/23

function add_patches_to_tuningdat(filetable_in, excel_rows)

excel_master_import;

wbar = waitbar(0,'Adding patch location info to tuning tables...');
for ifile = 1:nfiles
    tuning_file_name = [getfname(filetable.tuning_file{ifile}) '.mat'];
    cd(filetable.directory{ifile})
    load(tuning_file_name)
    tuningdat = get_roi_patch_location(tuningdat, filetable.patch_file{ifile}, filetable.reg_to_patches_file{ifile}); % add patch info
    tuningdat.input_files.patchdata_file = filetable.patch_file{ifile};
    tuningdat.input_files.patch_reg_file = filetable.reg_to_patches_file{ifile};
    delete(tuning_file_name) % delete so we have a backup in recycle bin
    save(tuning_file_name, 'tuningdat') % save tuning data with patch info
    clear tuningdat tuning_file_name
    try waitbar(ifile/nfiles,wbar); end
end
try close(wbar); end