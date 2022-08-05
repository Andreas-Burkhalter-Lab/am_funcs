%%%% BATCH_RENAME: replace 'string_to_remove' with 'string_to_insert' in the filenames of
%%%% all files in directory 'dirc'; copies of the originals will be put in
%%%% recycle bin
%   in:
%       1. directory
%       2. string to remove
%       3. string to insert
% last updated 2018/2/27

function batch_rename(dirc,string_to_remove,string_to_insert)

recycle('on');
if ~exist('dirc','var') || isempty(dirc)
    dirc = pwd;
end
cd(dirc);
filetable = struct2table(dir); % get all files in dirc

%get indices of files in the table to rename
namematch = find(~cellfun(@isempty,strfind(filetable.name,string_to_remove)));

files_to_rename = filetable.name(namematch);

for ind = 1:length(files_to_rename)
    [~,oldname,oldext] = fileparts(files_to_rename{ind});
    copyToDelete = [oldname '_copy' oldext];
    copyfile([dirc,filesep,oldname,oldext],[dirc,filesep,copyToDelete]); % make a copy
    delete(copyToDelete); % delete the copy
    newname = strrep(oldname,string_to_remove,string_to_insert);
    movefile([oldname oldext],[newname oldext]) % rename the file
end
