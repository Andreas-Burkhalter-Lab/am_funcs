%%% analyze rf mapping, tuning etc. for all sessions listed in master_xls, save results into multiple tuningdat files
% [filetable] = batch_full_session_analysis(filetable_in, excel_rows, savename, pars)
%       excel_rows are rows within the excel file, so the variable names row = 1, first actual data row = 2
%       savename = name to save each of the results files as; default is 'tuningdat'
%       pars = pars struct for get_tuning_curves
%%%%% updated 2019/1/30 on thermaltake

function [filetable] = batch_full_session_analysis(filetable_in, excel_rows, savename, pars)

excel_master_import;

savename = vardefault('savename','tuningdat'); 
pars = vardefault('pars',struct); 
pars.fit_tuning_functions = field_default(pars,'fit_tuning_functions',1);

wbar = waitbar(0,'Performing full session analysis...');
for ifile = 1:nfiles
    cd(filetable.directory{ifile})
    
    % set input files
    input_files.dff_file = filetable.dff_file{ifile};
    input_files.triggertiming_file = filetable.trigdata_file{ifile};
    input_files.regops_file = filetable.regops_file{ifile};
    input_files.scopetiming_file_rf = filetable.scopetiming_file_rf{ifile};
    input_files.scopetiming_file_dark = filetable.scopetiming_file_dark{ifile};
    input_files.scopetiming_file_tuning = filetable.scopetiming_file_tuning{ifile};
    input_files.stimdata_file_rf = filetable.stimdata_file_rf{ifile};
    input_files.stimdata_file_tuning = filetable.stimdata_file_tuning{ifile};
%     pars.pupil_data_file = filetable.pupil_data_file{ifile};
    input_files.pupil_data_file = [];
    input_files.patchdata_file = filetable.patch_file{ifile};
    input_files.patch_reg_file = filetable.reg_to_patches_file{ifile};
    
    % input defaults if empty
    input_files.pupil_data_file = field_default(input_files, 'pupil_data_file',[]);
    input_files.patchdata_file = field_default(input_files, 'patchdata_file',[]);
    input_files.patch_reg_file = field_default(input_files, 'patch_reg_file',[]);
    % run analysis
    tuningdat = full_session_analysis(input_files, pars);
    save(savename, 'tuningdat')
    clear input_files tuningdat
    try waitbar(ifile/nfiles,wbar); end
end
try close(wbar); end