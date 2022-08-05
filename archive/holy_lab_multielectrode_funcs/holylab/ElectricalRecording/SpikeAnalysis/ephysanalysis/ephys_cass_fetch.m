function ephyso = ephys_cass_fetch(ephysi,chanfile)
% EPHYS_CASS_FETCH: a temporary means for importing CASS sorting results
%
% NOTE: this is deprecated. Use the standard "ephysfetch" instead, it now
% has support for CASS results.  This is here only for supporting scripts
% which made use of it, new scripts should not use this.
%
% Syntax:
%   ephyso = ephys_cass_fetch(ephysi,chanfile)
% where
%   ephysi is your ephys structure array---you have to manually set the
%     "sortfile" field to the name of the directory (e.g., "sort2")
%     containing your sorting data.
%   chanfile is a structure array with the following fields:
%      channel: the channel number
%      file: the base name of the *.sorted file to use for the given
%        channel.
% and
%   ephyso is your output ephys structure with cell times loaded in.
%
% Note: because this is cheesy, you have to call it at the right
% time. You want to do this before you've called ephyssubrange.
% You also need to have the number of elements in your ephysi structure
% be equal to the number of files sorted by cass. There is no checking to
% make sure that these refer to the same files.
  
  nchans = length(chanfile);
  dirname = ephysi(1).sortfile;
  ephyso = ephysi;
  for i = 1:length(ephyso)
    ephyso(i).cellnums = [];
    ephyso(i).celltimes = {};
    ephyso(i).cellchandef = {};
  end
  for i = 1:nchans
    chanpath = [dirname filesep 'chan' num2str(chanfile(i).channel) ...
                filesep];
    if strcmp(chanpath(1),'/') % would mean a full path has been provided already
        eval(['load -mat ' chanpath chanfile(i).file '.sorted']);
    else
        eval(['load -mat ./' chanpath chanfile(i).file '.sorted']);
    end
    [ncells,nfiles] = size(chanclust);
    for j = 1:nfiles
      nextCell = 1;
      for k = 1:ncells
%        ephyso(j).cellnums(end+1) = chanfile(i).channel + k/100;
        ephyso(j).cellnums(end+1) = nextCell;
        nextCell = nextCell+1;
        ephyso(j).celltimes{end+1} = chanclust{k,j};
        ephyso(j).cellchandef{end+1} = chanfile(i).channel;
      end
      ephyso(j).sort_cass = 1;
    end
  end
  
      
  