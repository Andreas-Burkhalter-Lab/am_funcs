% 2008_06_02_nonparse_script.m

file_list = dirbyname('*.fmbx');
file_list_size = size(file_list,2);

if mod(file_list_size,2) == 0
    for idx = 1:file_list_size/2
        tic;
        nonsparse_fmb(file_list([idx idx+file_list_size/2]), 'glow')
        fprintf('File %s complete: ', file_list{idx});
        toc;
    end
else fprintf('Odd number of files, ''nonsparse_fmb'' must be used via command line in this folder');
    
end