%%% PrintFig.m
%  Loads and Prints Matlab figures listed in batch file.
%
%

batchfiledir = 'Z:\Data\Tempo\Batch Files\Binding\'
filename = input('Enter Batch File Name: ', 's');

status = ['Loading batch file: ' batchfiledir filename];
disp(status);

filename = [batchfiledir filename];
fid = fopen(filename);
output = 1;
print_fig = 0;
file_num = 1;

line = fgetl(fid);
while (line ~= -1)   
    if (line(1) ~= '%')
	    spaces = isspace(line);
		space_index = find(spaces);

		%get path / file
		PATH = line(1:space_index(1) - 1);
		FILE = line(space_index(1) + 1:space_index(2) - 1);
        FILE(1,end-2:end) = 'ccg'
        PATH = PATH(1,1:end-4);
        PATH = [PATH 'Analysis\Tuning\'];
        if (exist([PATH FILE]))  
            eval (['open ' PATH FILE]);    
            eval (['close ' PATH FILE]);    
        else
            sprintf('File %s not exist', FILE);    
        end %if ccg file exists       
    end % if (line(1) ~=...
  	line = fgetl(fid);
end %while...
fclose all;

