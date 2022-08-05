function snips = sortheader_readsnips(sh,index)
% SORTHEADER_READSNIPS: read spike waveforms
% Syntax:
%   snips = sortheader_readsnips(shc,index)
% where
%   shc is a channel-specific sortheader (see SORTHEADER_IMPORTCHAN);
%   index is a cell array of snippet indices, containing the list of
%     snippet waveforms you want to load (can be a vector instead of a
%     cell array if only one file is being supplied);
% and
%   snips is a cell array of length(shc), where each entry is a
%     sniplength-by-nsnips matrix of spike waveforms.  Snippets are
%     scaled to units of voltage, using scalemult and scaleoff.  (snips
%     will be a matrix instead of a cell array of matrices if index was
%     supplied as a vector.)
%
% See also: SORTHEADER_IMPORTCHAN, READINDEXSNIP.
  
% Copyright 2005 by Timothy E. Holy
  
  nfiles = length(sh);
  outcell = 1;
  if (nfiles == 1 && ~iscell(index))
    index = {index};
    outcell = 0;
  end
  for i = 1:nfiles
    if ~isfield(sh(i),'type') || strcmp(sh(i).type,'classic')
      fh = fopen(sh(i).fh);
      if (fseek(fh.fid,sh(i).snipsfpos,'bof') < 0)
	msg = ferror(fh.fid);
	error(['File ' fh.abspathsr fh.filename ' seek error: ' msg]);
      end
      snips{i} = readindexsnip(fh.fid,index{i},sh(i).sniprange, ...
			       sh(i).snipsprec);
      snips{i} = sh(i).scalemult*snips{i} + sh(i).scaleoff;
      fclose(fh);
    else
      switch sh(i).type
       case 'fit_component'
	error('Not implemented');
	%fullpath = [shin(i).dirname filesep 'chan' num2str(channel)];
	%timeampfilename = fullfile(fullpath,['time_amp' num2str(i)]);
        %load(timeampfilename);
      end
    end
  end
  if (~outcell & length(snips) == 1)
    snips = snips{1};
  end
