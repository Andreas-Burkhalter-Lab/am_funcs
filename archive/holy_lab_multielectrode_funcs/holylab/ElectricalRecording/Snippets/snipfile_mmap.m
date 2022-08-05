function mmap = snipfile_mmap(filename,channel)
% SNIPFILE_MMAP: read snippet file using memory mapping (for big files)
% Syntax:
%   mmap = snipfile_mmap(filename,channel)
% where
%   filename is a string containing the name of the file;
%   channel is an integer channel number you desire to read;
% and on output, mmap is a memmapfile object, with the following data
% fields:
%   mmap.data.time: contains the scan number
%   mmap.data.snip: contains the snippets
%   mmap.data.finetime: contains the fine times (may be absent, if
%     interp_times was false)
%   mmap.data.detpeak: contains the detection peak values
%
% See also: MEMMAPFILE.
  
  h = read_snip_header(filename);
  channelIndex = find(h.channels == channel);
  if isempty(channelIndex)
    error(['Channel ' num2str(channel) ' not in file ' filename]);
  end
  % Check to make sure times, snips, fine times, and detpeaks are
  % contiguous, as required by memmapfile
  % Construct a matrix of the following format:
  %   [fpos; expected_fpos; datasize]
  nsnips = h.numofsnips(channelIndex);
  snipsize = diff(h.sniprange)+1;
  ci2 = [channelIndex; channelIndex];
  % Times
  tfpos = h.timesfpos;
  tfpos = tfpos(:);
  offset = [tfpos(ci2); 4];  % [tfpos; tfpos; size(time)]
  % Snips
  offset = [offset, [h.snipsfpos(channelIndex); ...
                     offset(2,end)+nsnips*offset(3,end); ...
                     2*snipsize]];
  % Finetimes
  if ~isempty(h.finetimesfpos)
    offset = [offset, [h.finetimesfpos(channelIndex);...
                       offset(2,end)+nsnips*offset(3,end); ...
                       4]];
  end
  % Detpeaks
  offset = [offset, [h.detpeaksfpos(channelIndex); ...
                     offset(2,end)+nsnips*offset(3,end); ...
                     2]];
  if ~isequal(offset(1,:),offset(2,:))
    error(['Snippet data in ' filename ' are not contiguous, can''t use ' ...
                        'memory mapping'])
  end
  
  if ~isempty(h.finetimesfpos)
    mmap = memmapfile(filename,'offset',offset(1),'repeat',1,'format', {...
        'int32' [1 nsnips] 'time';...
        'int16' [snipsize nsnips] 'snip';...
        'single' [1 nsnips] 'finetime';...
        'int16' [1 nsnips] 'detpeak'});
  else
    mmap = memmapfile(filename,'offset',offset(1),'repeat',1,'format', {...
        'int32' [1 nsnips] 'time';...
        'int16' [snipsize nsnips] 'snip';...
        'int16' [1 nsnips] 'detpeak'});
  end