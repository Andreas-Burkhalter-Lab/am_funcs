function [set_data file_list] = gather_analysis(directory_list, varargin)
% gather_analysis loads analysis.mat files from the file list
% Syntax: set_data = gather_analysis(directory_list, varargin)
%
% Inputs: 
%         directory_list = a cell array of strings containing the
%                          directory(ies) where the desired analysis.mat
%                          files are located 
%         varargin: paired input arguments, currently limited to:
%                   >> fullpath: logical (default true)                                              
%              ** to do ** extract only certain subfields of analysis structures
%
% Outputs:
%    set_data:  a structure array of the structure type saved within 'analysis.mat'
%    file_list: optional output of the fullpath version of the files loaded 
%               (for cross-checking purposes)
% 
% See also SINGLE_ELECTRODE_EPHYS_ANALYZE, MULTI_ELECTRODE_EPHYS_ANALYZE, MEEA_APPLY

% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)
% Revision history:
%     2008_01_31: (JPM) wrote it
%     2009_08_13: cleaned help file, checked on function (JPM)

%% check for errors in file_list
if ~iscell(directory_list)
    if ~ischar(directory_list)
        errormsg = 'file_list must be a string or cell array of strings';
        error(errormsg);
    end
else
    for file_idx = 1:size(directory_list,2)
        if ~ischar(directory_list{file_idx})
            errormsg = 'all entries to file_list must be strings';
            error(errormsg);
        end
    end
end
%% check for varargin
if nargin > 1
    if ischar(varargin{1})
        switch varargin{1}
            case 'fullpath'
                if islogical(varargin{2})
                    fullpath = varargin{2};
                else
                    errormsg = 'fullpath must be a logical';
                    error(errormsg);
                end
            %%% add new option flags...
        end
    end
else
    fullpath = true;
end

file_list = directory_list;

%% if fullpath is FALSE
%     add current directory to all filenames to make the loading section
%     easy
if fullpath == false
    base = pwd;
    if iscell(file_list)
        for file_idx = 1:size(file_list,2)
            file_list(file_idx) = {[base filesep file_list{file_idx}]};
        end
    else
        file_list = {[base filesep file_list]};
    end
end
%% load up the analysis.mat files from each file in file_list
for file_idx = 1:size(file_list,2)
    if isempty(regexp(file_list{file_idx},'.mat', 'once'))
        if file_list{file_idx}(end) == '/'
            file_list{file_idx} = [file_list{file_idx} 'analysis.mat'];
        else
            file_list{file_idx} = [file_list{file_idx} filesep 'analysis.mat'];
        end
    end
    data = load(file_list{file_idx}, '-mat');
    for data_idx = 1:size(data.analysis,2)
        if ~exist('set_data', 'var')
            if iscell(data)
                set_data(1).analysis = {data.analysis(data_idx)};
            else
                set_data{1}.analysis = data.analysis(data_idx);
            end
        else
            if iscell(data)
                set_data(end+1).analysis = {data.analysis(data_idx)};
            else
                set_data{end+1}.analysis = data.analysis(data_idx);
            end
        end
    end

end