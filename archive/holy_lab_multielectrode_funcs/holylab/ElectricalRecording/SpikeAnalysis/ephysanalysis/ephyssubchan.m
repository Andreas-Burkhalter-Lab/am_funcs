function ephysout = ephyssubchan(ephysin,channels)
% EPHYSSUBCHAN: eliminate all but specified channels from ephys structure
% ephysout = ephyssubchan(ephysin,channels)
% 
% ephysin can be a structure array
%
% See also: EPHYS, EPHYSSUBRANGE, EPHYSSUBCELL.

  % Check that desired channels are present, and find channel indices
  nstruct = prod(size(ephysin));
  channels = channels(:)';
  nchannels = length(channels);
  chanIndex = cell(1,nstruct);
  for i = 1:nstruct
    [comchan,chanIndex{i}] = intersect(ephysin(i).channels,channels);
    if (length(comchan) < length(channels))
      error('Some desired channels are missing from input');
    end
  end
  chanIndex = cat(1,chanIndex{:});
  
  % Copy over any unmodified data
  fieldstomodify = {'channels','wave','envelope','snippets', ...
                    'sniptimes','snippeaks','snipthresh','celltimes',...
                    'cellnums','cellchandef'};
  ephysout = copyotherfields(ephysin,fieldstomodify);

  % Now deal with the fields that will be changed
  for i = 1:nstruct
    ephysout(i).channels = channels;
    if isfield(ephysin,'wave')
      ephysout(i).wave = ephysin(i).wave(chanIndex(i,:),:);
    end
    if isfield(ephysin,'envelope')
      echanIndex = sort([2*(chanIndex(i,:)-1)+1,2*chanIndex(i,:)]);
      ephysout(i).envelope = ephysin(i).envelope(echanIndex,:);
    end
    if isfield(ephysin,'sniptimes')
      ephysout(i).sniptimes = ephysin(i).sniptimes(chanIndex(i,:));
    end
    if isfield(ephysin,'snipindex')
      ephysout(i).snipindex = ephysin(i).snipindex(chanIndex(i,:));
    end
    if isfield(ephysin,'snippets')
      ephysout(i).snippets = ephysin(i).snippets(chanIndex(i,:));
    end
    if isfield(ephysin,'snippeaks')
      ephysout(i).snippeaks = ephysin(i).snippeaks(chanIndex(i,:));
    end
    if isfield(ephysin,'snipthresh')
      ephysout(i).snipthresh = ephysin(i).snipthresh(chanIndex(i,:));
    end
    if isfield(ephysin,'celltimes')
      % Find the cells which are defined entirely on the chosen channels
      keepcells = zeros(size(ephysin(i).cellchandef));
      for j = 1:length(keepcells)
        issel = ismember(ephysin(i).cellchandef{j},channels);
        keepcells(j) = all(issel);
      end
      keepcells = find(keepcells);
      % Copy data
      ephysout(i).cellnums = ephysin(i).cellnums(keepcells);
      ephysout(i).cellchandef = ephysin(i).cellchandef(keepcells);
      ephysout(i).celltimes = ephysin(i).celltimes(keepcells);
    end
    if (isfield(ephysin(i),'cellnums') & ischar(ephysin(i).cellnums))
      ephysout(i).cellnums = ephysin(i).cellnums;
    end
  end
