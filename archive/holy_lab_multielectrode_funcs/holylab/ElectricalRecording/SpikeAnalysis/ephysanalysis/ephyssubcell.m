function ephysout = ephyssubcell(ephysin,cellnums)
% EPHYSSUBCELL: eliminate all but specified cells from ephys structure
% ephysout = ephyssubcell(ephysin,cellnums)
% 
% ephysin can be a structure array
%
% See also: EPHYS, EPHYSSUBRANGE, EPHYSSUBCHAN.

  % Check that desired cell numbers are present, and find cellnum indices
  nstruct = numel(ephysin);
  cellnums = cellnums(:)';
  ncellnums = length(cellnums);
  cellIndex = cell(1,nstruct);
  for i = 1:nstruct
    [comcell,cellIndex{i}] = intersect(ephysin(i).cellnums,cellnums);
    if (length(comcell) < length(cellnums))
      error('Some desired cellnums are missing from input');
    end
  end
  % Determine which channels these cells showed up on, and compute the
  % channel index
  chanIndex = cell(1,nstruct);
  for i = 1:nstruct
    if (length(ephysin(i).cellchandef) > 1)
      chantmp = unique(cat(2,ephysin(i).cellchandef{:}));
    else
      chantmp = ephysin(i).cellchandef{1};
    end
    [comchan,chanIndex{i}] = intersect(ephysin(i).channels,chantmp);
    if (length(comchan) < length(chantmp))
      error('Some channels are missing for some selected cells');
    end
  end
  cellIndex = cat(1,cellIndex{:});
  chanIndex = cat(1,chanIndex{:});
  if isfield(ephysin(1),'sort_cass_chanfile')
    % Extract the cass channels
    cass_channels = unique(floor(cellnums));
    cass_chanIndex = findainb(cass_channels,[ephysin(1).sort_cass_chanfile.channel]);
  end
  
  % Copy over any unmodified data
  fieldstomodify = {'channels','wave','envelope','snippets', ...
                    'sniptimes','snippeaks','snipthresh','celltimes',...
                    'cellnums','cellchandef','sort_cass_chanfile','cellscanrange'};
  ephysout = copyotherfields(ephysin,fieldstomodify);

  % Now deal with the fields that will be changed
  for i = 1:nstruct
    ephysout(i).cellnums = cellnums;
    ephysout(i).channels = ephysin(i).channels(chanIndex(i,:));
    if isfield(ephysin,'wave')
      ephysout(i).wave = ephysin(i).wave(chanIndex(i,:),:);
    end
    if isfield(ephysin,'envelope')
      echanIndex = sort([2*(chanIndex(i,:)-1)+1,2*chanIndex(i,:)]);
      ephysout(i).envelope = ephysin(i).envelope(echanIndex,:);
    end
    ephysout(i).cellnums = ephysin(i).cellnums(cellIndex(i,:));
    ephysout(i).cellchandef = ephysin(i).cellchandef(cellIndex(i,:));
    if isfield(ephysin,'celltimes') && ~isempty(ephysin(i).celltimes)
      ephysout(i).celltimes = ephysin(i).celltimes(cellIndex(i,:));
    end
    if isfield(ephysin,'snipthresh')
      ephysout(i).snipthresh = ephysin(i).snipthresh(chanIndex(i,:));
    end
    if isfield(ephysin,'sniptimes')
      % The following code assumes cells are defined on a single
      % channel. No decisions have yet been made about how to handle the
      % representation in the more general case.
      %
      % First find out which snippets we need to keep. We have to
      % accumulate this across cells, as two cells may share a channel.
      keepindx = cell(1,length(ephysin(i).channels));
      for j = 1:ncellnums
        cItemp = find(ephysin(i).cellchandef{cellIndex(j)} == ephysin(i).channels);
        [ct,kItemp] = intersect(ephysin(i).sniptimes{cItemp},ephysin(i).celltimes{cellIndex(j)});
        keepindx{cItemp} = [keepindx{cItemp};kItemp(:)];
      end
      % Now (after we sort the index) we can throw away all but the
      % desired snippets
      for j = 1:length(keepindx)
        keepindx{j} = sort(keepindx{j});
        ephysout(i).sniptimes{j} = ephysin(i).sniptimes{j}(keepindx{j});
        if isfield(ephysin,'snipindex')
          ephysout(i).snipindex{j} = ...
              ephysin(i).snipindex{j}(keepindx{j});
        end
        if isfield(ephysin,'snippets')
          ephysout(i).snippets{j} = ...
              ephysin(i).snippets{j}(:,keepindx{j});
        end
        if isfield(ephysin,'snippeaks')
          ephysout(i).snippeaks{j} = ...
              ephysin(i).snippeaks{j}(keepindx{j});
        end
      end
    end
    if isfield(ephysin,'sort_cass_chanfile')
      ephysout(i).sort_cass_chanfile = ephysin(i).sort_cass_chanfile(cass_chanIndex);
    end
  end
