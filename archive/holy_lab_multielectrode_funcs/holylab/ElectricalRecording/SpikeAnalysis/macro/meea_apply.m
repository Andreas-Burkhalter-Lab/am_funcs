function passed = meea_apply(directory_list, varargin)
% meea_apply applies the most recent version of the seea_analysis to the directories in directory_list
% Syntax: passed = meea_apply(directory_list, varargin)
% 
% Required inputs:
% directory_list: a char array or cell array of strings containing
%                 directory names (full path)
% varargin: either (A) paired arguments for multi_electrode_ephys_analyze
%                  (B) a single struct containing meea option fields
% 
% Output:
%    passed: a true/false value that is true when the function has completed
%            all analyses successfully
%
% See also MULTI_ELECTRODE_EPHYS_ANALYZE

% Version History:
% 2008_06_23: Wrote it (JPM)
% 2009_08_13: cleaned up, checked (JPM)

%% Error checking
if ~iscell(directory_list)
    if ischar(directory_list)
        directory_list = {directory_list};
    else
        error('directory_list is not a char array or cell array of strings.');
    end
end

%% Execute script
home = pwd;

for f_idx = 1:size(directory_list, 2)
    cd(directory_list{f_idx});
    meea('directory', directory_list{f_idx}, varargin{1:end});
    fprintf('%s completed successfully.\n', directory_list{f_idx});
end

cd(home);
fprintf('%d files processed successfully.\n', size(directory_list,2));
passed = 1;

end