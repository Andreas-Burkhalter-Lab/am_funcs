function pathname = strip_filename(fullpathname)
%STRIP_PATH: returns only the path from an input string which includes directory info
%
% Copyright Julian P. Meeks 2008
    
fileseps = strfind(fullpathname, filesep);       % finds slashes by index within string
pathname = fullpathname(1:fileseps(end));
end
