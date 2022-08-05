%FSIZE: get size of file in KB, MB, or GB. 
%
%%% 1. 'filename' = name of file or directory to get size of
%%% 2. [optional] units = 'kb', 'mb' (default), or 'gb' 
%
% Output is filesize in specified units. If directory, get size of
% contents, including all levels of subdirectories.
%
%%% Last updated 2021/6/2 by Andrew Meier

function varargout = fsize(filename,units)

[a b c] = fileparts(filename);
if isempty(c)
    if isempty(a)
        d = struct2table(dir);
    else
        d = struct2table(dir(a));
    end
    match = strcmp(filename,cellfun(@getfname,d.name,'UniformOutput',false));
    if any(match)
        filename = fullfile(fileparts(filename),d.name{find(match,1)});
    end
end

if ~exist('units','var')
    units = 'mb';
end

switch units
    case 'kb'
        divideby = 1e3;
    case 'mb'
        divideby = 1e6;
    case 'gb'
        divideby = 1e9;
end

fhandle = dir(filename);
if isempty(fhandle)
    error('Did not find %s.',filename)
end

% following three commands copied from http://www.mathworks.com/matlabcentral/newsreader/view_thread/152712
if isdir(filename) % if it's a directory, get contents of all levels of subdirectories
    dirList = regexp(genpath(filename), pathsep, 'split');
    fileList = cellfun(@dir, dirList, 'Uniform', 0);
    filesize = sum(arrayfun(@(f) f.bytes, vertcat(fileList{:})));
else % if it's a file, not a directory
    filesize = fhandle.bytes;
end

sizeInSpecUnits = filesize/divideby; % convert to specified units

fprintf('%g%s = size of %s.\n',sizeInSpecUnits,units,filename)

if nargout>0
    varargout{1} = sizeInSpecUnits;
end