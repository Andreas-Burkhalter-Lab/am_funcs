function [tStatus, tOutput] = calculate_envelopes(pathstrs,basenames,merecfiles);

% calculate the envelope for a merec file
%
% SYNTAX: [tStatus, tOutput] = calculate_envelopes(pathstrs,basenames,merecfiles)
%           (any combination of these three inputs - including none -
%                                               can be provided in any order) 
%
% just a matlab front for the calenv utility so that it has a help and more
% flexible syntax; note that pathstrs can be a cell array of the same
% length as merecfiles or it can be a single string that works for all of
% them

%% deal with inputs: what do we actually have?
current_directory = pwd;
getmerecs = 1;
getbasenames = 1;
getpathstrs = 1;
if nargin == 1
    if (iscell(pathstrs) && isstrfind(pathstrs{1},'.merec')) | isstrfind(pathstrs,'.merec')
        merecfiles = pathstrs;
        getmerecs = 0;
    elseif ~((iscell(pathstrs) && isstrfind(pathstrs{1},filesep)) | isstrfind(pathstrs,filesep))
        basenames = pathstrs;
        getbasenames = 0;
    else
        getpathstrs = 0;
    end
elseif nargin == 2
    temp = basenames;
    % figure out what the first input is...
    if (iscell(pathstrs) && isstrfind(pathstrs{1},'.merec')) | isstrfind(pathstrs,'.merec')
        merecfiles = pathstrs;
        getmerecs = 0;
    elseif ~((iscell(pathstrs) && isstrfind(pathstrs{1},filesep)) | isstrfind(pathstrs,filesep))
        basenames = pathstrs;
        getbasenames = 0;
    else
        getpathstrs = 0;
    end
    % figure out what the second input is...
    if (iscell(temp) && isstrfind(temp{1},'.merec')) | isstrfind(temp,'.merec')
        merecfiles = temp;
        getmerecs = 0;
    elseif ~((iscell(temp) && isstrfind(temp{1},filesep)) | isstrfind(temp,filesep))
        basenames = temp;
        getbasenames = 0;
    else
        pathstrs = temp;
        getpathstrs = 0;
    end
elseif nargin == 3
    getmerecs = 0;
    getbasenames = 0;
    getpathstrs = 0;
end

%% deal with inputs: fill in the gaps and reformat

% get merecfiles and pathstrs
if getmerecs && getpathstrs
    merecfiles = dirbytime('*.merec');
    nFiles = length(merecfiles);
    pathstrs = cell(1,nFiles);
    [pathstrs{:}] = deal(pwd);
elseif getmerecs
    temp = pathstrs;
    pathstrs = cell(1,0);
    merecfiles = cell(1,0);
    if ~iscell(temp)
        temp = {temp};
    end
    nFilesSoFar = 0;
    for nthPath = 1:length(temp);
        cd(temp{nthPath})
        merectemp = dirbytime('*.merec');
        [merecfiles{(nFilesSoFar+1):(nFilesSoFar+length(merectemp))}] = deal(merectemp{:});
        [pathstrs{(nFilesSoFar+1):(nFilesSoFar+length(merectemp))}] = deal(temp{nthPath});
        nFilesSoFar = nFilesSoFar+length(merectemp);
    end
elseif getpathstrs
    if iscell(merecfiles)
        merecfiles = {merecfiles};
    end
    pathstrs = cell(1,length(merecfiles));
    [pathstrs{:}] = deal(pwd);
end
% get basenames
if getbasenames
    basenames = cell(1,length(merecfiles)); % fill in later rather than put in loop now
end

%% do it
for nthFile = 1:length(merecfiles)
    if ~strcmp(pathstrs{nthFile}(end),filesep)
        pathstrs{nthFile} = [pathstrs{nthFile} filesep];
    end
    if getbasenames
        basenames{nthFile} = merecfiles{nthFile}(1:(end-6));
    end
    [tStatus, tOutput]=system(['calenv -s 100 -o ' pathstrs{nthFile} basenames{nthFile} ...
                    '.env '  merecfiles{nthFile}]);
end
                
                