function shout = sortheader_importchan_AM(shin,channel,file_override)
% SORTHEADER_IMPORTCHAN: focus on a specific channel for sorting
% This function also loads in the spike time and peak amplitude
%

%%% AM updated 10/20/15 on vivid: add option to specify new name/location of
%%% .ssnp file (in case it was moved/copied from original location)

% Syntax:
%   shc = sortheader_importchan(sh,channel)
% where
%   sh is a structure array of sortheaders (see SNIPFILE2SORTHEADER);
%   channel is the channel number;
% and
%   shc is the channel-specific sortheader, which contains the additional
%     fields "sniptimes" (a vector of spike times, in scan #s) and
%     "peakHeight" (measured in the units of spike waveforms, and is to
%     be compared relative to "thresh").
%
% See also: SNIPFILE2SORTHEADER.
  
% Copyright 2005 by Timothy E. Holy
  
  shout = shin;
  nfiles = length(shin);
  for i = 1:nfiles
    chanindx = find(shin(i).channels == channel);
    if isempty(chanindx)
      error(['Channel ' num2str(channel) ' not included in file ' ...
             shin(i).fh.abspathstr shin(i).fh.filename]);
    end
    shout(i).channels = channel;
    shout(i).thresh = shin(i).thresh(:,chanindx);
    if ~isfield(shin,'type') || strcmp(shin(i).type,'classic')
      fn = {'numofsnips','timesfpos','snipsfpos','finetimesfpos','detpeaksfpos'};
      for j = 1:length(fn)
        if ~isempty(shin(i).(fn{j}))
          shout(i).(fn{j}) = shin(i).(fn{j})(chanindx);
        end
      end
      
      % Get sniptimes (assume we can get these all in memory)
      if shin(i).fh.uselfs
        error('Not implemented')
      else
        [fh,msg] = fopen(shin(i).fh);
        
       %% AM added 7/25/15
        if exist('file_override','var') && ~isempty(file_override)
            [temppath tempname tempext] = fileparts(file_override);
            fh.abspathstr = temppath;
            if ~isempty(fh.abspathstr)
                fh.abspathstr = [fh.abspathstr filesep]; % restore filesep removed by fileparts
            end
            fh.filename = [tempname tempext];
            fh.fid = fopen(file_override);
            shout.fh = fh; % shout inherits the new overriding filename
        end
        %%
        
        if (fh.fid < 0)
          error(['Error opening ' fh.abspathstr fh.filename ': ' msg]);
        end
        if (fseek(fh.fid,shout(i).timesfpos,'bof') < 0)
          msg = ferror(fh.fid);
          error(['File ' fh.abspathstr fh.filename ' seek error: ' msg]);
        end
        shout(i).sniptimes = fread(fh.fid,shout(i).numofsnips,...
          ['*' shout(i).timesprec]);
        if isempty(shout(i).sniptimes)
          shout(i).sniptimes = zeros(0,1);
        end
        if (fseek(fh.fid,shout(i).detpeaksfpos,'bof') < 0)
          msg = ferror(fh.fid);
          error(['File ' fh.abspathstr fh.filename ' seek error: ' msg]);
        end
        shout(i).peakHeight = fread(fh.fid,shout(i).numofsnips,...
          ['*' shout(i).detpeaksprec]);
        fclose(fh);
      end
    else
      switch shin(i).type
        case 'fit_component'
          shout(i).numofsnips = shin(i).numofsnips(chanindx);
          fullpath = [shin(i).dirname filesep 'chan' num2str(channel)];
          timeampfilename = fullfile(fullpath,['time_amp' num2str(i)]);
          tt = load(timeampfilename,'fc_spiketimes');
          shout(i).sniptimes = tt.fc_spiketimes(:);
          % For peakHeight, import the min/max info
          if (i == 1)
            projfilename = fullfile(fullpath,'proj.mat');
            tt = load(projfilename,'snipmm');
            snipmm = tt.snipmm;
          end
          thismm = snipmm{i};
          [mxmm,mmindx] = max(abs(thismm));
          mmindx = sub2ind([2 length(mmindx)],mmindx,1:length(mmindx));
          shout(i).peakHeight = thismm(mmindx)';
          % Also import the rank
          tt = load(timeampfilename,'fc_ranks');
          if isfield(tt,'fc_ranks')
            shout(i).fc_rank = tt.fc_ranks;
          end
      end
    end
  end
  