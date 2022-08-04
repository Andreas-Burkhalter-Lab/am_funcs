%%%% import master excel file listing recording sessions to analyze and associated files
%   -master_xls and excel_rows should be input arguments from the parent function
%   -if filetable_in is inputted as a table rather than a filename, output that table without adjustment
%
% updated 2020/3/25 on thermaltake

if ~istable(filetable_in) % if filetable_in is a filename rather than an already-processed table variable
    filetable = excel2table(filetable_in); % load master file list
    excel_rows = vardefault('excel_rows',2:height(filetable)+1);
    matlab_table_rows = excel_rows - 1; % compensate for variable names row taking up one place in excel file
    filetable = filetable(matlab_table_rows,:); % only keep rows specified by excel_rows
    filetable.analyze_plane(isnan(filetable.analyze_plane)) = 0; 
    filetable = filetable(logical(filetable.analyze_plane),:); % keep only planes that have been marked for analysis
    filetable.analyze_plane = []; 
    nfiles = height(filetable);

    % make into cells for compatibility
    if isnumeric(filetable.scopetiming_file_rf) && all(isnan(filetable.scopetiming_file_rf))
        filetable.scopetiming_file_rf = cell(height(filetable),1);
    end
    if isnumeric(filetable.stimdata_file_rf) &&  all(isnan(filetable.stimdata_file_rf))
        filetable.stimdata_file_rf = cell(height(filetable),1);
    end
    if isnumeric(filetable.scopetiming_file_dark) && all(isnan(filetable.scopetiming_file_dark))
        filetable.scopetiming_file_dark = cell(height(filetable),1);
    end
elseif istable(filetable_in)
    filetable = filetable_in;
end