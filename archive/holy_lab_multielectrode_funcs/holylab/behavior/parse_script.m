% 2008_06_02_parse_script.m

file_list = dirbyname('*.fmb');

for idx = 1:size(file_list, 2)
    tic;
    parse_fmb(file_list{idx})
    fprintf('File %s complete: ', file_list{idx});
    toc;
end
