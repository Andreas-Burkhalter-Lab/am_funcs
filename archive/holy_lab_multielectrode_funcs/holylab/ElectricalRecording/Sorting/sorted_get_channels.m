function channels = sorted_get_channels(dirname,ext)
% SORTED_GET_CHANNELS: get the channel #s for sorted data
% Syntax:
%   channels = sorted_get_channels(dirname)
%   channels = sorted_get_channels(dirname,ext)
%
% Scans the directory "dirname" for subdirectories with the name "chanX"
% where X is a channel #.  It returns the list of channel numbers of
% directories with extension ".mat".  Alternatively, you can supply a
% different extension yourself.
%
% Doing things this way, rather than relying on the channel list in
% options_autosort, insures that users can copy over just the
% "interesting" directories and expect things to still work properly.
%
% See also: PREPARE_CHANFILE.

% Copyright 2005 by Timothy E. Holy
  
% Something to think about: what to do if a cell is defined using
% multiple channels? This could probably be encoded in the directory
% name, e.g. chan41_53
  
  if (nargin < 2)
    ext = '.mat';
  end

  chandirs = dir([dirname filesep 'chan*']);
  chandirIndex = [chandirs.isdir];  % keep only the directories
  chandirnames = {chandirs(chandirIndex).name};
  channels = [];
  for i = 1:length(chandirnames)
    if (length(chandirnames{i}) > 4)
      tchannel = str2num(chandirnames{i}(5:end));
      if ~isempty(tchannel)
        % Directory has the name 'chanX', where X is a channel #
        % Is the directory empty, or does it have files in it?
        thisChanDir = dir([dirname filesep chandirnames{i} filesep ...
                            '*' ext]);
        if ~isempty(thisChanDir)
          channels(end+1) = tchannel;
        end
      end
    end
  end
  channels = sort(channels);
