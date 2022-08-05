function filename = strip_path(fullpathname)
%STRIP_PATH: returns only the filename from an input string which includes directory info
%
% Copyright Julian P. Meeks 2007
    
fileseps = strfind(fullpathname, filesep);       % finds slashes by index within string
filename = fullpathname(fileseps(end)+1:length(fullpathname));
end
