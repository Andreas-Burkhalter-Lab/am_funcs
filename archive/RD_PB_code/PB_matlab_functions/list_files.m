function fname = list_files(data_directory,filter);
% This function lists files with the 'filter'  in 'data_directory'
% Inputs: data_directory in strings, e.g. '/home/mukesh'
%       : filter in strings, e.g. '*.cnt'

finfo = dir(strcat(data_directory,'/',filter));
fname = char(finfo.name);
disp([num2str(size(fname,1)) ' files']);
