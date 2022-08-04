function [varargout] = sorthead_from_raw_all( varargin )
%SORTHEAD_FROM_RAW_ALL Saves snips and sniptimes from raw data extracted by autosort
%for all .merec files in the working directory. 
%  Run this function after running autosort_all.
%  Assumes that directories exist with the same names as all .merec files
%  in the working directory; each of these directories should contain the
%  'overview' file created by autosort. These directories and files are
%  created by autosort_all. 

filenames = dir('*.ssnp');
[placeholder filenames] = cellfun(@fileparts,extractfield(filenames,'name'),'UniformOutput',false);

for i = 1:length(filenames)
    filenames{i}        
    clear sorthead
    load(strcat(filenames{i},'/overview.mat'))
    snipstruct = sorthead_from_raw(sorthead);
    save(strcat(filenames{i},'/snipstruct.mat'),'snipstruct')
end


end

