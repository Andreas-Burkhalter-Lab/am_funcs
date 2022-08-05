function ephysout = ephysfetch(ephysin,fldnames)
% EPHYSFETCH: load electrophysiology data
% Syntax:
%   ephysout = ephysfetch(ephysin,fieldnames)
%
% ephysin and ephysout are ephys structures; see EPHYS.
% fieldnames is a cell array of the names of ephys fields, from among the
%   following choices:
%     'stimulus','wave','envelope','snippets','sniptimes',
%     'snippeaks','cellnums','celltimes'
%
% Alternatively, you can get information that is not an ephys field:
%     numsnips = ephysfetch(ephysin,'numsnips')
% where numsnips is a matrix, where numsnips(i,j) is the number of
% snippets on the jth channel in the _whole snippet file_ corresponding
% to the ith ephys structure.
% You can't mix&match this syntax with the field-type nomenclature.
%     
%
% The names of the data files must already be set (e.g., 'stimulusfile'
% must be set if you are loading in information about 'stimulus').
%
% If you don't want to specify the cell numbers manually, you may set the
% 'cellnums' field to 
%  'all': to get all cells, irrespective of the 'channels' field
%  'allonchan': to get all cells, as long as they were defined on
%     channels included in the channels field.
% Fetching 'cellnums' will replace these strings with the real cell numbers.
% (In progress: 'channels' gets an 'all' option)
%
% See also: EPHYSTOFETCH, EPHYSSUBRANGE, EPHYSSUBCHAN, EPHYSSUBCELL, EPHYS.

% Note: in progress: generalizing to input structure arrays

  if ischar(fldnames)
    fldnames = {fldnames};
  end
  nfields = length(fldnames);

  % Check special syntaxes
  if (nfields == 1)
    switch fldnames{1}
     case 'numsnips'
      options = struct;
      if isfield(ephysin,'snipfilemachfmt')
        options.machfmt = ephysin(1).snipfilemachfmt;
        options.headertype = 'snip';
      end
      for i = 1:length(ephysin)
        [h,fid] = readheader(ephysin(i).snipfile,options);
        chanindx = findainb(ephysin(i).channels,h.channels);
        ephysout(i,:) = h.numofsnips(chanindx);
      end
      return
    end
  end
  
  % Regular syntax
  allfields = {'stimulus','wave','envelope','snippets', ...
               'sniptimes','snippeaks','cellnums','celltimes'};
  % Check that we know where to get all the data
  for i = 1:nfields
    if ~any(strmatch(fldnames{i},allfields,'exact'))
      error(['Don''t know about field ',fldnames{i}]);
    end
    % Several types of data come from snippet files,
    % so need a special check here
    if strncmp(fldnames{i},'snip',4)
      if ~isfield(ephysin,'snipfile')
        error('Don''t know where to get snippet information');
      end
    elseif (strncmp(fldnames{i},'cell',4))
      if ~isfield(ephysin,'sortfile')
        error(['Don''t know source for loading information for field ', ...
               fldnames{i}]);
      end
    elseif ~isfield(ephysin,[fldnames{i},'file'])
      error(['Don''t know source for loading information for field ', ...
             fldnames{i}]);
    end
  end
  
  % Copy any old data
  ephysout = ephysin;
  
  % stimulus
  i = strmatch('stimulus',fldnames,'exact');
  if ~isempty(i)
    % Figure out which stimulus files we need to read
    stimulusfiles = unique({ephysout.stimulusfile});
    for j = 1:length(stimulusfiles)
      eindx = strmatch(stimulusfiles{j},{ephysout.stimulusfile},'exact');
      [tempfiles,tempstim] = readvlv(stimulusfiles{j});
      for k = 1:length(eindx)
        cur = eindx(k);
        i = strmatch(ephysout(cur).basefilename,tempfiles,'exact');
        if (length(i) ~= 1)
          error('File names do not match uniquely, or at all');
        end
        % Truncate the range
        lindx = [1,find(tempstim{i}(2,:) <= ephysout(cur).scanrange(1))];
        rindx = find(tempstim{i}(2,:) <= ephysout(cur).scanrange(2));
        rng = [lindx(end):rindx(end),rindx(end)];
        ephysout(cur).stimulus = tempstim{i}(:,rng);
        ephysout(cur).stimulus(2,[1 end]) = ephysout(cur).scanrange;
      end
    end % stimulusfiles
  end % stimulus
  
  % waveform
  i = strmatch('wave',fldnames,'exact');
  if ~isempty(i)
    for i = 1:length(ephysout)
      options = struct;
      if isfield(ephysout(i),'wavefilemachfmt')
        % Not necessary except for legacy files
        options.machfmt = ephysout(i).wavefilemachfmt;
        options.headertype = 'ai';
      end
      [h,fid] = readheader(ephysout(i).wavefile,options);
      if (ischar(ephysout(i).channels) & strcmp(ephysout(i).channels, ...
                                                'all'))
        comchan = h.channels;
        chanIndex = 1:length(comchan);
      else
        [comchan,chanIndex] = intersect(h.channels,ephysout(i).channels);
      end
      if (length(comchan) < length(ephysout(i).channels))
        error('Some desired channels are missing from file');
      end
      fseek(fid,(ephysout(i).scanrange(1)-1)*h.numch*2,'cof');
      nscans = diff(ephysout(i).scanrange) + 1;
      prec = 'int16';
      if isfield(ephysout(i),'wavefileprec')
        prec = ephysout(i).wavefileprec;
      end
      temp = fread(fid,[h.numch nscans],prec);
      fclose(fid);
      ephysout(i).wave = temp(chanIndex,:)*h.scalemult + h.scaleoff;
    end
  end
  
  % envelope
  i = strmatch('envelope',fldnames,'exact');
  if ~isempty(i)
    for i = 1:length(ephysout)
      options = struct;
      if isfield(ephysout(i),'envelopefilemachfmt')
        options.machfmt = ephysout(i).envelopefilemachfmt;
        options.headertype = 'env';
      end
      [h,fid] = readheader(ephysout(i).envelopefile,options);
      if (ischar(ephysout(i).channels) & strcmp(ephysout(i).channels, ...
                                                'all'))
        chanIndex = 1:length(h.channels);
      else
        %[comchan,chanIndex] = intersect(h.channels,ephysout(i).channels);
        chanIndex = findainb(ephysout(i).channels,h.channels);
      end
      %chanIndex = sort([2*chanIndex-1,2*chanIndex]); % min/max pairs
      chanIndex = [2*chanIndex-1; 2*chanIndex];
      chanIndex = chanIndex(:)';
      fseek(fid,round((ephysout(i).scanrange(1)-1)/h.decimate)*h.numch*4,'cof');
      nscans = diff(ephysout(i).scanrange) + 1;
      ndec = round(nscans/h.decimate);
      prec = 'int16';
      if isfield(ephysout(i),'envelopefileprec')
        prec = ephysout(i).envelopefileprec;
      end
      temp = fread(fid,[2*h.numch ndec],prec);
      fclose(fid);
      ephysout(i).envelope = temp(chanIndex,:)*h.scalemult + h.scaleoff;
      ephysout(i).envelopedecimate = h.decimate;
    end
  end
  
  % spiketimes, snippets, and spikepeaks
  tnames = {'sniptimes','snippets','snippeaks'};
  [comnames,indx] = intersect(strvcat(tnames),strvcat(fldnames), ...
                              'rows');
  if ~isempty(indx)
    %if any(indx == 3)  % snippeaks
    %  error('Not implemented')
    %end
    % Determine which files need to be read, and process
    % them in parallel for efficiency
    ufiles = unique({ephysout.snipfile});
    for k = 1:length(ufiles)
      sindx = strmatch(ufiles{k},{ephysout.snipfile},'exact');
      options = struct;
      if isfield(ephysout,'snipfilemachfmt')
        options.machfmt = ephysout(sindx(1)).snipfilemachfmt;
        options.headertype = 'snip';
      end
      [h,fid] = readheader(ufiles{k},options);
      allchannels = unique([ephysout(sindx).channels]);
      [comchan,chanIndex] = intersect(h.channels,allchannels);
      if (length(comchan) < length(allchannels))
        % Check w/ user that the missing channels aren't ones that
        % should have no snippets on them, e.g. the stim. channel
        missingchans = setdiff(allchannels,h.channels);
        buttonname = questdlg(['Channels ',num2str(missingchans),...
                    ' are missing from snippet file. Proceed anyway?'],'','Yes','No','Yes');
        if strcmp(buttonname,'No')
          error('Some desired channels are missing from file');
        end
      end
      if (prod(size(h.thresh)) > length(h.channels))
        [ephysout(sindx).snipthresh] = deal(h.thresh(:,chanIndex) * h.scalemult + h.scaleoff);
      else
        % Legacy files
        [ephysout(sindx).snipthresh] = deal(h.thresh(chanIndex) * h.scalemult + h.scaleoff);
      end
      [ephysout(sindx).snippolarity] = deal(h.options.polarity);
      intervals = cat(1,ephysout(sindx).scanrange);
      % Load times & compute index numbers, if necessary
      % (we need to load the times, even to load just the snippets,
      %  for indexing purposes)
      loadtimes = 0;
      if (any(indx == 1 | indx == 2) & ~isfield(ephysout,'snipindex'))
        loadtimes = 1;
        sindxsub = sindx;
      elseif any(indx == 1 | indx == 2)
        % Check to see if snipindex exists but is empty
        needsloading = zeros(1,length(sindx));
        loadindx = ephystofetch(ephysout(sindx),'snipindex');
        needsloading(loadindx) = 1;
        % Check to see if all-empty cell array
        nlindx = setdiff(1:length(sindx),loadindx);
        noload = sindx(nlindx);
        needsloadingcheck = ones(1,length(noload));
        for j = 1:length(noload)
          for jj = 1:length(ephysout(noload(j)).channels)
            if ~isempty(ephysout(noload(j)).snipindex{jj})
              needsloadingcheck(j) = 0;
            end
          end
        end
        if any(needsloadingcheck)
%          keyboard    % commenting out and putting warning in instead:
           warning('There are intervals with no sniptimes....')
        end
        needsloading(nlindx) = needsloadingcheck;
        if any(needsloading)
          loadtimes = 1;
          sindxsub = sindx(find(needsloading));
        end
      end
      if loadtimes
        % OK, let's load the times
        for i = 1:length(chanIndex)
          fseek(fid,h.timesfpos(chanIndex(i)),'bof');
          ttemp = fread(fid,h.numofsnips(chanIndex(i)),'int32');
          itemp = timesininterval(ttemp,intervals);
          for j = 1:length(sindxsub)
            echanIndex = find(ephysout(sindxsub(j)).channels == comchan(i));
            if ~isempty(echanIndex)
              ephysout(sindxsub(j)).sniptimes{echanIndex} = ttemp(itemp{j});
              ephysout(sindxsub(j)).snipindex{echanIndex} = itemp{j};
            end
          end
        end
      elseif any(indx == 1)   % snipindex is already defined, so use it
        for i = 1:length(chanIndex)
          fseek(fid,h.timesfpos(chanIndex(i)),'bof');
          ttemp = fread(fid,h.numofsnips(chanIndex(i)),'int32');
          for j = 1:length(sindx)
            echanIndex = find(ephysout(sindx(j)).channels == comchan(i));
            if ~isempty(echanIndex)
              ephysout(sindx(j)).sniptimes{echanIndex} = ...
                    ttemp(ephysout(sindx(j)).snipindex{echanIndex});
            end
          end
        end
      end
      % Load snippets
      if any(indx == 2)
        [ephysout(sindx).sniprange] = deal(h.sniprange);
        prec = 'int16';
        if isfield(ephysout,'snipfileprec')
          prec = ephysout(sindx(1)).snipfileprec;
        end
        for i = 1:length(chanIndex)
          for j = 1:length(sindx)
            echanIndex = find(ephysout(sindx(j)).channels == comchan(i));
            if ~isempty(echanIndex)
              fseek(fid,h.snipsfpos(chanIndex(i)),'bof');
              ephysout(sindx(j)).snippets{echanIndex} = h.scaleoff + h.scalemult * ...
                readindexsnip(fid,ephysout(sindx(j)).snipindex{echanIndex},h.sniprange,prec);
            end
          end
        end
      end % Done with snippets
      % Load snippeaks
      if any(indx == 3)
        for i = 1:length(chanIndex)
          fseek(fid,h.detpeaksfpos(chanIndex(i)),'bof');
          ptemp = fread(fid,h.numofsnips(chanIndex(i)),'int16');
          for j = 1:length(sindx)
            echanIndex = find(ephysout(sindx(j)).channels == comchan(i));
            if ~isempty(echanIndex)
              ephysout(sindx(j)).snippeaks{echanIndex} = h.scaleoff + h.scalemult * ...
                ptemp(ephysout(sindx(j)).snipindex{echanIndex});
            end
          end
        end
      end % done with snippeaks  
      fclose(fid);
    end % Done with current snipfile
  end % Done with loop over different snipfiles
  
  % cells: cellnums and celltimes
  if (~isempty(strmatch('cellnums',fldnames,'exact')) | ...
      ~isempty(strmatch('celltimes',fldnames,'exact')))
    loading_celltimes = ~isempty(strmatch('celltimes',fldnames,'exact'));
    if(isfield(ephysin, 'sort_method') && isequal(ephysin(1).sort_method, 'cluster nav'))

       % calc the "global" dominant channel for each cell
       % NOTE: may not be accurate when there're cycles having diff sides
       for sideIndex=1:length(ephysout(1).sortfile)
          navFile=ephysout(1).sortfile{sideIndex};
          nav=load(navFile, '-mat');
          preFile=replace_filename(navFile, nav.preClusterNavFile);
          pre=load(preFile,'fitFiles', '-mat');
          allChannels=[]; 
          clear fitcomp
          for fitFileIndex=1:length(pre.fitFiles)
             fitFile=replace_filename(navFile, pre.fitFiles{fitFileIndex});
             fitcomp{fitFileIndex}=load(fitFile, 'peakChan', 'componentFile', '-mat');
             if(isempty(allChannels))
                templatefile=replace_filename(navFile, fitcomp{fitFileIndex}.componentFile);
                template=load(templatefile, 'channels', '-mat');
                allChannels=template.channels;
             end
          end
          fitcomp=[fitcomp{:}];
          peakChan=[fitcomp.peakChan];
          
          clusters=agglabel(nav.spikeClust+1);
          
          for cellIndex=1:length(clusters)
             peakChanCur=peakChan(clusters{cellIndex});
             if(isempty(peakChanCur))
                peakChannels{sideIndex}(cellIndex)=NaN;
             else
                [peakChanCurUni, tt, indices]=unique(peakChanCur);
                indices_cluster=agglabel(indices);
                [tt, i]=max(cellfun(@length, indices_cluster));
                peakChannels{sideIndex}(cellIndex)=allChannels(peakChanCurUni(i(1)));
             end
          end % for, each cell
       end %for, each side

       
       for ephysIndex=1:length(ephysout)
          ephysout(ephysIndex).celltimes={};
          ephysout(ephysIndex).cellnums=[];
          ephysout(ephysIndex).cellchandef={};
          
          for sideIndex=1:length(ephysout(ephysIndex).sortfile)
             % put side loop inside ephys loop is slow but allow more
             % flexible config such as discarding some field/side for
             % specific cycles
             navFile=ephysout(ephysIndex).sortfile{sideIndex};
             nav=load(navFile, '-mat');
             preFile=replace_filename(navFile, nav.preClusterNavFile);
             pre=load(preFile,'fitFiles','cumspikes','-mat');
             fitFile=[ephysout(ephysIndex).basefilename '.fitcomp'];
             fitFileIndex=strmatch(fitFile, pre.fitFiles, 'exact');
             spikeClust=nav.spikeClust(pre.cumspikes(fitFileIndex)+1:pre.cumspikes(fitFileIndex+1));
             if isempty(spikeClust)
               continue
             end
             fitcomp=load(replace_filename(navFile, fitFile), 'spiketimes', '-mat');
             if(~isequal(nav.clusterIDs{1}, 'noise'))
                error('i''m lost --- there must be a noise cluster');
             end
             spike_labels=spikeClust+1; % due to that agglabel() requires labels>=1
             clusters=agglabel(spike_labels);
             celltimes=clusters(nav.locked==1);
             cellNames=nav.clusterIDs(nav.locked==1);
             cellchandef=peakChannels{sideIndex}(nav.locked==1);
             cellnums=NaN(1, length(cellNames));
             for cellIndex=1:length(cellNames)
                cellnums(cellIndex)=str2num(cellNames{cellIndex}(2:end))+sideIndex/1000; % drop the leading c
                celltimes{cellIndex}=fitcomp.spiketimes(celltimes{cellIndex}); % map col indices to scan numbers
             end
             
             ephysout(ephysIndex).celltimes=[ephysout(ephysIndex).celltimes celltimes];
             ephysout(ephysIndex).cellnums=[ephysout(ephysIndex).cellnums cellnums];
             ephysout(ephysIndex).cellchandef=[ephysout(ephysIndex).cellchandef num2cell(cellchandef)];
          end % for, each side (i.e. field of the array)
       end % for, each ephys struct
    elseif isfield(ephysin,'sort_cass') && ephysin(1).sort_cass
      % Load from CASS-sorted files
      dirname = unique({ephysin.sortfile});
      if (length(dirname) > 1)
        error('Must use a consistent sorting directory!');
      end
      dirname = dirname{1};
      % Load the overview.mat file; check to see whether the user tried to
      % load data from template-sorted files
      ovview = load([dirname filesep 'overview']);
      if isfield(ovview,'isFake') && ovview.isFake
        error('This does not work for template-based sorting; use load_fake_channel_celltimes instead');
      end
      % Verify that user specified the sortfile to use
      if ~isfield(ephysin,'sort_cass_chanfile')
        error('You must supply the list of channels and the sorting files to use.');
      end
      % Check that the chanfile is the same for all elements of the structure
      for ephysIndex = 2:length(ephysin)
        if ~isequal(ephysin(ephysIndex).sort_cass_chanfile,ephysin(1).sort_cass_chanfile)
          error('The chanfile has to be the same for all elements in the structure');
        end
      end
      chanfile = ephysin(1).sort_cass_chanfile;
      % Match ephys referents to their sort results
      sorted_snipfiles = {};
      for fileIndex = 1:length(ovview.sorthead)
        sorted_snipfiles{end+1} = ovview.sorthead(fileIndex).fh.filename;
      end
      sorting_file_index = findainb({ephysin.snipfile},sorted_snipfiles);
      % Wipe out any pre-existing cell info
      for epIndex = 1:length(ephysout)
        ephysout(epIndex).cellnums = [];
        ephysout(epIndex).celltimes = {};
        ephysout(epIndex).cellchandef = {};
        ephysout(epIndex).cellscanrange = {};
      end
      % Now get to work!
      for chanfileIndex = 1:length(chanfile)
        this_cass_chan = chanfile(chanfileIndex).channel;
        chanpath = [dirname filesep 'chan' num2str(this_cass_chan) filesep];
        % Process the sort_info to get the number of cells on this channel
        load([chanpath chanfile(chanfileIndex).file]);
        n_cells = length(unique(sort_info.landmarkClust)) - 1;
        if (n_cells == 0)
          continue;   % No cells on this channel
        end
        % Add the new cell#s to the list
        % Note: the reason cells are numbered channel.number (i.e., 3.02
        % for the second cell on channel 3) is that it's the only way of
        % insuring that cell numbers do not change if the considered
        % channels change. While you might think that you can make a unique
        % identifier by combining cellchandef and the cellnumber,
        % cellchandef can be a vector (useful for multielectrode sorting)
        % and is meant to indicate the real "electrical" channels rather
        % than any "fake channel" used in sorting. Thus, because the "fake
        % channel" is potentially useful, we need to encode it somehow. One
        % could add a separate field to ephys (casschannel), but then one
        % would have to start checking multiple fields for identity,
        % ephysplotparams would need a new field, etc. So it seems easiest
        % to do it this way.
        cellListIndex = length(ephysout(1).cellnums) + (1:n_cells);
        cellnums = this_cass_chan + (1:n_cells)/100;
        for ephysIndex = 1:length(ephysout)
          ephysout(ephysIndex).cellnums(cellListIndex) = cellnums;
        end
        % Process the sort_info to get breakpoints
        file_intervals = timemarkers2intervals(sort_info.timeMarker,...
          [ovview.sorthead.nscans],n_cells);
        % Load the spike times
        if loading_celltimes
          channel_sorting_results = load('-mat',[chanpath chanfile(chanfileIndex).file '.sorted']);
        end
        for cellIndex = 1:n_cells
          for ephysIndex = 1:length(ephysin)
            this_file = sorting_file_index(ephysIndex);
            % Intersect the valid intervals & scanranges
            this_intervals = IntersectIntervals(ephysin(ephysIndex).scanrange,...
              file_intervals{cellIndex,this_file});
            this_cell_index = cellListIndex(cellIndex);
            ephysout(ephysIndex).cellscanrange{this_cell_index} = ...
              this_intervals;
            if loading_celltimes
              % Constrain the celltimes to be within these intervals
              t = channel_sorting_results.chanclust{cellIndex,this_file};
              if (isequal(this_intervals,ephysin(ephysIndex).scanrange) ...
                  && ~isempty(t) && t(1) >= this_intervals(1) && t(end) < this_intervals(2))
                % Can keep all of them
                ephysout(ephysIndex).celltimes{this_cell_index} = t;
              else
                timesindex = timesininterval(t,this_intervals);
                timesindex = sort(cat(1,timesindex{:}));
                ephysout(ephysIndex).celltimes{this_cell_index} = t(timesindex);
                %             indxkeep = false(size(t));
                %             for intervalIndex = 1:size(this_intervals,1)
                %               indxkeep = indxkeep | (t >= this_intervals(intervalIndex,1) & t < this_intervals(intervalIndex,2));
                %             end
                %             ephysout(ephysIndex).celltimes{this_cell_index} = t(indxkeep);
              end
            end
            % Record the chandef
            ephysout(ephysIndex).cellchandef{this_cell_index} = this_cass_chan;
          end
        end
      end
    else
      % Find the different sorting files
      sortfiles = unique({ephysin.sortfile});
      % Process each sorting file
      for i = 1:length(sortfiles)
        indx = strmatch(sortfiles{i},{ephysin.sortfile},'exact');
        try
          load(sortfiles{i})
        catch
          err = lasterror;
          fprintf('There was an error loading a sorting file called %s;\nmight you have forgotten to set sort_cass to true?\n',sortfiles{i});
          rethrow(err);
        end
        % Generate the basename from the .ssnp filenames
        for j = 1:length(spikefiles)
          [pathstr,name] = fileparts(spikefiles{j});
          sortspikefiles{j} = name;
        end
        % Concatenate data to a more convenient format:
        % catchans(cellindex,fileindex)
        catchans = cat(1,chanclust{:});
        % Determine the cell numbers
        lastnum = 0;
        for k = 1:length(chanclust)
          ncells = size(chanclust{k},1);
          cellnumtemp{k} = lastnum + [1:ncells]';
          channumtemp{k} = repmat(channels(k),ncells,1);
          lastnum = lastnum + ncells;
        end
        filecellnums = cat(1,cellnumtemp{:});
        channelsbycell = cat(1,channumtemp{:})';
        % For each element of ephysout, find the corresponding
        % file in the list of sorting files, get elements of
        % chanclust, and copy them over
        for j = 1:length(indx)
          % Find the right file
          indxs = strmatch(ephysout(indx(j)).basefilename,sortspikefiles,'exact');
          if (length(indxs) ~= 1)
            error('Error matching sort files and basefilename');
          end
          % Compare channels, if necessary
          cellnums = ephysout(indx(j)).cellnums;
          if ischar(cellnums)
            if strcmp(cellnums,'allonchan')
              ikeep = ismember(channelsbycell,ephysout(indx(j)).channels);
              cellnums = find(ikeep);
            elseif strcmp(cellnums,'all')
              cellnums = 1:length(channelsbycell);
            else
              error('Unrecognized option to cellnums');
            end
          end
          ephysout(indx(j)).cellnums = cellnums;
          if loading_celltimes
            % Copy over the data, restricting the spike times to within the
            % scan range (fix TEH 2004-06-01)
            % (Note on the face of it this seems inefficient, since
            % timesininterval can be more efficient if you ask for several
            % intervals from one set of times. Look into changing this if
            % performance becomes an issue.)
            for k = 1:length(cellnums)
              ephysout(indx(j)).cellchandef{k} = channelsbycell(cellnums(k));
              keepindx = timesininterval(catchans{cellnums(k),indxs},ephysout(indx(j)).scanrange);
              ephysout(indx(j)).celltimes{k} = catchans{cellnums(k),indxs}(keepindx{1});
            end
            %ephysout(indx(j)).celltimes = catchans(cellnums,indxs);
          end % if celltimes as well as cellnums
        end % loop over elements using same sorting file
      end % loop over sorting files
    end % non-CASS
  end % celltimes
