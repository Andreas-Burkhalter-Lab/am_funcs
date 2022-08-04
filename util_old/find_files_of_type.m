 %%%% find all files with a given extension with working directory; dot not required
 % matching_files = find_files_of_type(extension_string)
 
 %%%% updated 2018/12/22 on thermaltake
 
 function matching_files = find_files_of_type(extension_string)
 
 dd = struct2table(dir);
 dd = dd(3:end,:); %% first 2 entries are . and ..
 extlength = length(extension_string);
 deleterow = cellfun(@(x)length(x),dd{:,1}) < extlength; %% delete short filenames
    dd(deleterow,:) = [];
 extension_match = cell2mat(cellfun(@(x)strcmp(x(end-[extlength-1]:end),extension_string), dd{:,1}, 'UniformOutput',0)); 
 matching_files = dd.name(extension_match);