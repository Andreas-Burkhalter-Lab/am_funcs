function strHeader = snipfile_append_channel(fid,strHeader,channel,varargin)
% SNIPFILE_APPEND_CHANNEL: append a channel to a snippet file
% Syntax:
% To write an entire set of snippets:
%   strHeader = snipfile_append_channel(fid, strHeader, channel, times,
%                                       snips, finetimes, detpeaks)
% To write a subset of snippets from a channel too big to hold in memory
% at once:
%   strHeader = snipfile_append_channel(fid, strHeader, channel, mmap,
%                                       keepflag)
% With this syntax, only times(keepflag), snips(:,keepflag), etc are
% written to disk.  keepflag is a logical vector.  This syntax is useful
% when memory mapping, see below.
%
% Note: if snips is too big to fit in memory at once, use SNIPFILE_MMAP.
%
% See also: SNIPFILE_MMAP.

  % Parse input args
  use_mmap = 0;
  if strcmp(class(varargin{1}),'memmapfile')
    use_mmap = 1;
    s = varargin{1};
    nsnips = length(s.data.time);
    if (length(varargin) > 1)
      keepflag = varargin{2};
    else
      keepflag = logical(ones(1,nsnips));
    end
  else
    % reorganize input to make it look like mmap
    s.data.time = varargin{1};
    s.data.snip = varargin{2};
    s.data.finetime = varargin{3};
    s.data.detpeak = varargin{4};
    nsnips = length(s.data.time);
  end

  % Just check to make sure this is a valid channel
  channels = str2num(key2value(strHeader,'channel list'));
  channelIndex = find(channels == channel);
  if isempty(channelIndex)
    error('Channel not present in file');
  end

  fseek(fid,0,'eof');   % Go to the end of the file
  fieldname = {'time','snip','finetime','detpeak'};
  fieldformat = {'int32','int16','float','int16'};
  fpos = [0 0 0 0];   % We have to save the file position for each type
                      % of data
  for i = 1:length(fieldname)
    % Write the data
    if (isfield(s.data,fieldname{i}) && ~isempty(s.data.(fieldname{i})))
      fpos(i) = ftell(fid);
      if use_mmap
        % Write in blocks, to prevent overfilling memory
        blocksize = 5000;
        nblocks = ceil(nsnips/blocksize);
        for k = 1:nblocks
          startIndex = (k-1)*blocksize+1;
          endIndex = min(nsnips,startIndex+blocksize-1);
          writeIndex = find(keepflag(startIndex:endIndex))+(startIndex-1);
          fwrite(fid,s.data.(fieldname{i})(:,writeIndex),fieldformat{i});
        end
      else
        % We have it in memory, so just write directly
        fwrite(fid,s.data.(fieldname{i}),fieldformat{i});
      end
    end
  end

  % Fix up the header
  headerfieldnames = {'timesfpos','snipsfpos','finetimesfpos', ...
                      'detpeaksfpos'};
  for i = 1:length(headerfieldnames)
    value = str2num(key2value(strHeader,headerfieldnames{i}));
    value(channelIndex) = fpos(i);
    strHeader = update_value(strHeader, headerfieldnames{i}, ...
                             num2str(value));
  end
  numofsnips = str2num(key2value(strHeader,'numofsnips'));
  num_real_snips = nsnips;
  if use_mmap
    num_real_snips = sum(keepflag);
  end
  numofsnips(channelIndex) = num_real_snips;
  strHeader = update_value(strHeader, 'numofsnips', num2str(numofsnips));
  
  % Write the header to disk
  update_header(fid,strHeader);
  
    
  