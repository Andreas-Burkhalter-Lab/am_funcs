function chanIndex = CheckAIFilesConsistent(filenames,channels,optionsin)
% CHECKAIFILESCONSISTENT: verify that all AnalogInput files contain desired info
% chanIndex = CheckAIFilesConsistent(filenames,channels,options)
% where
%    filenames is a cell array of file names
%    channels (optional) is a vector of channel #s that must have been
%      recorded in all files
%    options is a set of options to check for greater consistency:
%      options(1) = 1: require scan rates to be identical (default 0)
%      options(2) = 1: require that all scalemults be identical (default 0)
%  Fill these in as far as you need to go, and defaults will be provided
%  for the rest.
%
% The chanIndex output is a matrix, each column of which is the
% the channel index for the given file (i.e. recchannels(chanIndex) = channels)
if (nargin < 3)
        optionsin = [];
end
if (nargin < 2)
        channels = [];
end
nfiles = length(filenames);
if (nfiles == 0)
        return
end
% Fill in defaults on any unspecified options
default_options = [0 0];
options = default_options;
options(1:length(optionsin)) = optionsin(1:length(optionsin));

if (nargout > 0)
        chanIndex = zeros(length(channels),nfiles);
end
% First read in all the relevant header data
filechans = cell(1,nfiles);
scanrate = zeros(1,nfiles);
scalemult = zeros(1,nfiles);
for i = 1:nfiles
        h = ReadAIHeader(filenames{i});
        filechans{i} = h.channels;
        scanrate(i) = h.scanrate;
        scalemult(i) = h.scalemult;
end
% Now check that all requested channels were recorded
for i = 1:nfiles
        if (any(setdiff(channels,filechans{i})))
                error('Not all requested channels were recorded');
        end
        % If desired, compute the channel indices
        if (nargout > 0)
                [c,ia,ib] = intersect(channels,filechans{i});
                chanIndex(:,i) = ib(ia)';
        end
end
% Check all other parameters
% scanrate:
if (options(1))
        if (any(scanrate ~= scanrate(1)))
                error('Not all scanrates are the same');
        end
end
% scalemult:
if (options(2))
        if (any(scalemult ~= scalemult(1)))
                error('Not all scalemults are the same');
        end
end
