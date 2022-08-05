function passed = seea_apply(directory_list, seea_func)
% seea_mix_apply applies the most recent version of the seea_analysis to the directories in directory_list
% Syntax: passed = seea_apply(directory_list, seea_func)
% 
% Required inputs:
% directory_list: a char array or cell array of strings containing
%                 directory names (full path)
% seea_func: a seea_function supplied as a function handle (e.g. @seea_mix)
% 
% Output:
% passed: a true/false value reaching true when the list is complete
%

% Version History:
% 2008_06_23: Wrote it (JPM)

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
    seea('directory', directory_list{f_idx}, 'call_func', seea_func);
end

cd(home);
fprintf('%d files processed successfully.\n', size(directory_list,2));
passed = 1;

end