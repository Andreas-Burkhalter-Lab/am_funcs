function snipminmax = sortheader_snipminmax(sh,index)
% SORTHEADER_SNIPMINMAX: read the min/max of each snippet
% Syntax:
%   snipminmax = sortheader_snipminmax(shc)
%   snipminmax = sortheader_snipminmax(shc,index)
% where
%   shc is a channel-specific sortheader (see SORTHEADER_IMPORTCHAN);
%   index (optional) is a cell array of snippet indices, containing the
%     list of snippet waveforms you want to load (can be a vector instead
%     of a cell array if only one file is being supplied). If omitted, it
%     defaults to processing all snippets.
% and
%   snipminmax is a cell array of length(shc), where each entry is a
%     2-by-nsnips matrix, where the first row contains the minium value
%     of each spike waveform, and the 2nd row contains the maximum value
%     of each spike waveform.  Snippets are scaled to units of voltage,
%     using scalemult and scaleoff.  (snipminmax will be a matrix instead
%     of a cell array of matrices if index was supplied as a vector.)
%
% See also: SORTHEADER_IMPORTCHAN, SORTHEADER_READSNIP, READINDEXSNIP.
  
% Copyright 2006 by Timothy E. Holy
  
  nfiles = length(sh);
  outcell = 1;
  if (nargin < 2)
    for i = 1:nfiles
      index{i} = 1:length(sh(i).sniptimes);
    end
  end
  if (nfiles == 1 && ~iscell(index))
    index = {index};
    outcell = 0;
  end
  % Process snippets in chunks, since we might ask for more minmax
  % values than we have room for in complete snippets
  blocksize = 5000;  % process in blocks of 5000 snippets
  for i = 1:nfiles
    fh = fopen(sh(i).fh);
    if (fseek(fh.fid,sh(i).snipsfpos,'bof') < 0)
      msg = ferror(fh.fid);
      error(['File ' fh.abspathsr fh.filename ' seek error: ' msg]);
    end
    offset = 1;
    snipmmblock = {};
    while (offset <= length(index{i}))
      blockend = min(offset-1+blocksize,length(index{i}));
      sniptmp = readindexsnip(fh.fid,index{i}(offset:blockend),...
			      sh(i).sniprange,sh(i).snipsprec);
      snipmmtmp = min(sniptmp);
      snipmmtmp(2,:) = max(sniptmp);
      snipmmblock{end+1} = snipmmtmp;
      offset = offset+blocksize;
      % readindexsnip assumes you're at the beginning of this channel's
      % snippets, so we have to go back
      if (fseek(fh.fid,sh(i).snipsfpos,'bof') < 0)
        msg = ferror(fh.fid);
        error(['File ' fh.abspathsr fh.filename ' seek error: ' msg]);
      end
    end
    snipminmax{i} = sh(i).scalemult*cat(2,snipmmblock{:}) + sh(i).scaleoff;
    fclose(fh);
  end
  if (~outcell & length(snipminmax) == 1)
    snipminmax = snipminmax{1};
  end
