% this version analyzes all runs in the directories specified rather than
% just a specified subset of runs

% get .htb filenames from all subdirectories of the working directory into a list for getTempoTrialData
%%% last edited 4/27/16 on msi

% topdir should contain folders containing .htb files 
n_subjects = 14; % must match number of folders (subjects) within topdir to be analyzed)
topdir = 'C:\Users\AM\Documents\Lab\recordings'; 
dirstruct = dir(topdir);
folderlist = {dirstruct.name};
folderlist = folderlist(3:n_subjects+2); % first two directors are always . and ..

% get .htb filenames from each individual folder
allfiles = cell(1,1);
folder_prefixes = cell(1,1);
for ind = 1:length(folderlist)
    thisfolder = folderlist{ind};
    fileshere = dir(thisfolder);
    fileshere = {fileshere(3:end).name}; % first two directors are always . and ..
    ishtb = strfind(fileshere,'.htb');
    ishtb = ~cellfun(@isempty,ishtb);
    htbshere = fileshere(ishtb);
    allfiles = [allfiles; htbshere'];
    folder_prefixes = [folder_prefixes; repmat({thisfolder},length(htbshere),1)];
end

allfiles = allfiles(2:end);
folder_prefixes = folder_prefixes(2:end);

fseps = repmat({filesep},length(allfiles),1);
files_with_dirs = strcat(folder_prefixes,fseps,allfiles);

% below section eliminates loom receding runs because they cannot yet be
% read by looming speed tuning analysis script
loom_recede_listfile = 'C:\Users\AM\Documents\Lab\recordings\loom_recede_files.xlsx';
[junk loom_recede_files] = xlsread(loom_recede_listfile);
loom_recedes_to_delete = [];
for ind = 1:length(allfiles)
    if any(strcmp(getfname(allfiles{ind}),loom_recede_files))
        loom_recedes_to_delete = [loom_recedes_to_delete; ind];
    end
end
allfiles(loom_recedes_to_delete) = [];
folder_prefixes(loom_recedes_to_delete) = [];
fseps(loom_recedes_to_delete) = [];
files_with_dirs(loom_recedes_to_delete) = [];























