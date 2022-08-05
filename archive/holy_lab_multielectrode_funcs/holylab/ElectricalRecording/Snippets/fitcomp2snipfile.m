function fitcomp2snipfile(dirs,basenames)
% FITCOMP2SNIPFILE: write a channel-based ssnp file from fitcomp
% Syntax:
%   fitcomp2snipfile(dirs,basenames)
% where
%   dirs is a cell array of directory names (e.g.: {'left','right'}), or a
%     single directory name if you only want 1 (e.g., 'left');
%   basenames is a cell array of base file names to process
%     (e.g.: {'cycle1','cycle2','cycle3'}).
% You should run this in the directory containing the .merec files.
% The .fitcomp files have to be in the directories listed in dirs.
%
% This function will write .ssnp files that can be imported using
% ephysfetch.
%
% See also: FIT_COMPONENTS, EPHYSFETCH.

% Copyright 2009 by Timothy E. Holy

  if ~iscell(dirs)
    dirs = {dirs};
  end
  if ~iscell(basenames)
    basenames = {basenames};
  end
  n_dirs = length(dirs);
  n_files = length(basenames);
  for dirIndex = n_dirs:-1:1
    tmpl(dirIndex) = load('-mat',[dirs{dirIndex} filesep dirs{dirIndex} ...
      '.templates']); %#ok<AGROW>
  end
  for fileIndex = n_files:-1:1
    hdr(fileIndex) = readheader([basenames{fileIndex} '.merec']); %#ok<AGROW>
  end
  % Get the channel list
  channels = [tmpl.channels];
  n_chans = length(channels);
  % Compute the thresholds in A/D units
  for dirIndex = n_dirs:-1:1
    medv = (tmpl(dirIndex).medv - hdr(1).scaleoff)/hdr(1).scalemult;
    thrsh = tmpl(dirIndex).thresh/hdr(1).scalemult;
    tmp = [medv-thrsh medv+thrsh]';
    thresh{dirIndex} = tmp(:)'; %#ok<AGROW>
  end
  thresh = int16(cat(2,thresh{:}));
  
  default_sniphdrc = {'[snippet]',...
		      'source=fitcomp',...
          ['thresh=' sprintf('%g ',thresh)],...
		      sprintf('snipbeginoffset=%d',-n_chans+1),...
		      sprintf('snipendoffset=%d',n_chans),...
		      'interptimes=0',...
		      'interpsnips=0',...
          'numofsnips=',...
          'timesfpos=',...
          'snipsfpos=',...
          'finetimesfpos=',...
          'detpeaksfpos='};
  default_sniphdr = sprintf('%s\n',default_sniphdrc{:});

  for fileIndex = 1:n_files
    fprintf('Working on %s...\n',basenames{fileIndex});
    channel = cell(1,n_dirs);
    times = cell(1,n_dirs);
    snips = cell(1,n_dirs);
    detpeaks = cell(1,n_dirs);
    for dirIndex = 1:n_dirs
      fc = load('-mat',[dirs{dirIndex} filesep basenames{fileIndex} ...
        '.fitcomp']);
      channel{dirIndex} = tmpl(dirIndex).channels(fc.peakChan);
      times{dirIndex} = fc.spiketimes;
      sniptmp = permute(fc.minmax,[3 2 1]);
      sniptmp = sniptmp(:,:);
      snips{dirIndex} = sniptmp';  % indexed by time, chan-then-minmax
      detpeaks{dirIndex} = (fc.peakVal - hdr(fileIndex).scaleoff)/hdr(fileIndex).scalemult;
    end
    channel = cat(2,channel{:});
    times = cat(2,times{:});
    snips = cat(2,snips{:});
    detpeaks = cat(2,detpeaks{:});
    % Convert snips to A/D units
    snips = (snips - hdr(fileIndex).scaleoff)/hdr(fileIndex).scalemult;

    % Create the header. We strip the MEREC and replace it with SNIPPET
    % (for the header type). We also append the default header. Finally, we
    % add file-specific information.
    strHdr = sprintf('SNIPPET\n%s\n%s',hdr(fileIndex).wholeheader(7:end),default_sniphdr);
    strHdr = [strHdr sprintf('\nsnippet input file=%s\n',[basenames{fileIndex} '.merec'])];
    strHdr = update_value(strHdr,'channel list',num2str(channels));
    outname = [basenames{fileIndex} '.ssnp'];
    [fid,msg] = fopen(outname,'w');
    if (fid < 1)
      error(['Working on file ' basenames{fileIndex} ': ' msg]);
    end
    count = fwrite(fid,strHdr,'char');
    if (count < length(strHdr))
      error(['Error writing ' outname]);
    end
    headersize = str2double(key2value(strHdr,'header size'));
    padsize = headersize - ftell(fid);
    fwrite(fid,zeros(1,padsize),'uint8');
    if (ftell(fid) ~= headersize)
      error(['Error padding header for file ' outname]);
    end
    for chanIndex = 1:n_chans
      thisChannel = channels(chanIndex);
      keepFlag = channel == thisChannel;
      strHdr = snipfile_append_channel(fid,strHdr,thisChannel,...
        times(keepFlag),...
        snips(:,keepFlag),...
        zeros(1,sum(keepFlag)),...
        detpeaks(keepFlag));
    end
    fclose(fid);
  end
				       
    
    
      