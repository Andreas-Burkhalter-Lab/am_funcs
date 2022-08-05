function fit_components(componentFile, fileToFit, optionsForFit)
% FIT_COMPONENTS: identify spikes & resolve overlaps in multielectrode data
% Syntax:
%   fit_components(componentFile, fileToFit, options)
% where
%   componentFile: a string containing the name of a *.fine_cluster or
%     *.fine_cluster_svd file, containing either templates or 
%     components of templates (see SVDTEMPLATES)
%   fileToFit: a string containing the name of a .merec file
%   options: may have the following fields:
%     saveDirectory: if specified, output files are saved to this
%       directory
%     fileToSave: the filename used to save the fitting results (spike
%       times and component amplitudes). By default, this will be
%       generated from the fileToFit by replacing the extension with
%       ".fitcomp".  If this file is specified, it overrides the
%       saveDirectory information, so be sure to supply full path
%       details if needed.
%     isSaveRes (default true): if true, a _residual.merec file is
%       created, containing the waveform of the residual after fitting.
%     residualFile: the filename used to save the residual. By default,
%       this will be generated from the fileToFit by replacing the
%       extension with "_residual.merec". As with fileToSave,
%       saveDirectory will be ignored in defining the path if you supply
%       this field, so be sure to include any necessary path info.
%     max_components (default Inf): use only up to the given number of
%       waveform components in the fitting. This option is activated only if
%       the componentFile has a variable "lambda," indicating it's the
%       result of an SVD.
%     blockSize (default 2^14): the number of scans to process at
%       once. It's best to use a power of 2, to facilitate the fourier
%       transforms.  Be aware that a block of this size is used for each
%       (multichannel) component, so if you have a lot of components you
%       can quickly eat up a lot of memory.  Edges are handled gracefully
%       by including overlap, so this parameter should have no impact on
%       the final results.
%     isMergeAdjacent (default true): if true, take triggers that happen
%       shortly after a previous trigger and ignore it. This helps insure
%       that cells with supra-threshold spike waveforms on multiple
%       channels trigger only a single event.
%     adjacentNScans (default 5): the maximum number of scans between
%       triggers for them to be considered "adjacent".  Events separated
%       by more scans than this are considered separate events.
%     fitIterMax (default 2): if larger than 1, the residual after
%       fitting will be examined, and areas that are above threshold will
%       be checked to see if new events need to be inserted.  Note that this
%       way even events closer than adjacentNScans can be created.  This
%       parameter controls the number of times this process occurs; for
%       fitIterMax = 2, two fits (but only the second one refined by
%       examination of the residual for additional events) will be
%       performed.  The results of the final fit are what gets recorded.
%
% See also: COMPONENTS_FROM_WAVEFORM.
  
% Copyright 2007-2008 by Timothy E. Holy
   
  %% Initialization
  if(nargin==2)
    optionsForFit=struct;
  end
  optionsForFit=fillOptions(optionsForFit);
  fu = fit_utilities;

  templateVars=load(componentFile, '-mat');
  components=templateVars.templates;
  if isfield(templateVars,'lambda')
    max_components = min(optionsForFit.max_components,size(components,3));
    components = components(:,:,1:max_components);
  end
  channels=templateVars.channels;
  thresh=templateVars.thresh;
  medv=templateVars.medv;
  sniprange = templateVars.sniprange;

  if isfield(templateVars,'rawClusterFilename')
    rawClusterFile=templateVars.rawClusterFilename;
    tt=load(rawClusterFile, 'options', '-mat');
    rawClusterOptions=tt.options;
  else
    rawClusterOptions.polarity = templateVars.polarity;
  end
  if isfield(optionsForFit,'debug')
    rawClusterOptions.debug = optionsForFit.debug;
  end

  nChannels = length(channels);

  % file to save results
  if(isfield(optionsForFit, 'fileToSave'))
    fileToSave=optionsForFit.fileToSave;
  elseif isfield(optionsForFit,'saveDirectory')
    [pathstr,main_name] = fileparts(fileToFit);
    fileToSave = fullfile(optionsForFit.saveDirectory,[main_name '.fitcomp']);
  else
    fileToSave=replace_extension(fileToFit, '.fitcomp');
  end

  isSaveRes=optionsForFit.isSaveRes;
  if(isSaveRes)
    if isfield(optionsForFit,'residualFile')
      resFile = optionsForFit.residualFile;
    elseif isfield(optionsForFit,'saveDirectory')
      [pathstr,main_name] = fileparts(fileToFit);
      resFile = fullfile(optionsForFit.saveDirectory,[main_name '_residual.merec']);
    else
      resFile=replace_extension(fileToSave, '_residual.merec');
    end
    % We'll also use a key to make sure that the residual file is being
    % paired with the right set of times&amplitudes
    rand('state', sum(100*clock));  % generate _different_ keys on each matlab startup
    residual_key = round(1e8*rand(1,1));
  end

  %% Set up blocksize, i.e., the size of chunk loaded into memory
  blockSize=optionsForFit.blockSize;
  % We need to have a little bit of overlap between blocks, so that we
  % don't miss spikes on the boundary
  snipLen = diff(sniprange)+1;
  advanceSize = blockSize - 2*snipLen; % put buffer of snipLen on each edge

  %% Pre-calculations for the fit
  if (size(components,1) ~= nChannels)
    error('Size of components does not match channel list');
  end
  fprintf('\nProcessing file %s\n',fileToFit);
  fprintf('Calculating component overlaps...');
  c_overlap = component_overlaps(components);
  fprintf('done\n');
  components = double(components);

  %% Prepare the waveform loading & residual-saving
  memm = merecmm(fileToFit,'tovolts',true,'contiguous',true);
  header=memm.header;
  if isSaveRes
    fidRes = ftcmp_openresidual(resFile,header,channels,residual_key);
    % Calculate things needed to convert between floats and int16
    merecinfo.minsample = header.minSample;
    merecinfo.slope = (header.voltageMax-header.voltageMin)/ ...
      (header.maxSample-header.minSample);
    merecinfo.offset = header.voltageMax - ...
      merecinfo.slope*header.maxSample;
  end
  nBlocks=ceil(header.nscans/advanceSize);
  nBlocks = min(nBlocks,optionsForFit.maxBlocks);  % useful for debugging (terminate early)
  all_amp = cell(1,nBlocks);
  all_time = cell(1,nBlocks);
  all_peak = cell(1,nBlocks);
  all_chan = cell(1,nBlocks);
  all_pp = cell(1,nBlocks);
  all_err = cell(1,nBlocks);

  hfig = -1;  % for graphical debugging

  %%
  % Go through the raw file in blocks, extracting best-fit parameters for
  % each block. We'll save the results as we go (writing a residual.merec
  % file during each block, the fitting parameters at the end as a .mat
  % file). The edges of each block are contaminated, but by using
  % overlapping blocks we insure that each scan is in the "interior" of
  % exactly one block.
  %
  for idxBlock=1:nBlocks
    %% Load & preprocess the data
    fprintf('.');
    if mod(idxBlock,10) == 0
      fprintf('%d%% done\n',round(100*idxBlock/nBlocks));
    end
    scanNumFrom=(idxBlock-1)*advanceSize+1;
    scanNumTo  =min(scanNumFrom+blockSize-1, header.nscans);
    wave=memm(channels, [scanNumFrom scanNumTo]);
    [d,N] = size(wave);
    wave=wave-repmat(medv,1,N); % remove any bias diff. from 0

    %% Find the events (based on peak voltage)
    [eventTime,peakVal,peakChan] = findpeaks_multichan(wave,thresh,rawClusterOptions);
    % Kill events that are on the edge of the time window
    eventKeep = fu.flag_times_in_range(eventTime,[1 N],sniprange);
    eventTime = eventTime(eventKeep);
    peakVal = peakVal(eventKeep);
    peakChan = peakChan(eventKeep);

    % We have to do the fitting twice; once using the given times, and
    % then once with some new times added back in for those events that, in
    % the residual, are still too big
    for fitIter = 1:optionsForFit.fitIterMax
      %% Fit the amplitudes
      eventTimeShift = eventTime + sniprange(1);
      b = wdecomp_amplitude_rhs(wave,components,eventTimeShift);
      % Solve for the component amplitudes
      a = fu.solve_for_amplitudes(eventTimeShift,b,c_overlap);
      % Calculate the residual
      wave_res = fu.calculate_residual(wave,eventTimeShift,components,a);
      if (fitIter < optionsForFit.fitIterMax)
        %% Try to improve the fit by inserting additional events
        % Areas with big error can correspond to missed spikes---these
        % can occur due to overlaps
        [newTimes,newPV,newPC] = findpeaks_multichan(wave_res,thresh,rawClusterOptions);
        eventKeep = newTimes+sniprange(1) > 0 & ...
          newTimes+sniprange(2) <= N;
        [eventTime,uIndex] = unique([eventTime newTimes(eventKeep)]);
        peakVal = [peakVal newPV(eventKeep)];
        peakVal = peakVal(uIndex);
        peakChan = [peakChan newPC(eventKeep)];
        peakChan = peakChan(uIndex);
      end
    end
    %% Compute the min/max and error associated with each event
    [pp,err] = fu.minmax_using_residual(wave_res,eventTimeShift,components,a);

    %% Graphical debugging
    if optionsForFit.debug
      if ~ishandle(hfig)
        hfig = figure;
      else
        figure(hfig);
      end
      clf(hfig);
      haxdebug = subplot(2,1,1);
      plot(scanNumFrom:scanNumTo,wave');
      title('Original')
      haxdebug(2) = subplot(2,1,2);
      plot(scanNumFrom:scanNumTo,wave_res');
      title('Residual')
      ylim = get(haxdebug,'YLim');
      ylim = cat(1,ylim{:});
      ylim = [min(ylim(:,1)) max(ylim(:,2))];
      set(haxdebug,'YLim',ylim);
      slidwincmenu(haxdebug([2 1]));
      pause
      c = get(gcf,'CurrentCharacter');
      if strcmp(c,'q')
        % Turn off all debugging
        optionsForFit.debug = false;
        rawClusterOptions.debug = false;
        delete(hfig)
        drawnow
      elseif strcmp(c,'Q')
        % Just turn of findpeaks debugging
        rawClusterOptions.debug = false;
      end
    end

    %% Preserve these items for the end when we save the results
    % But keep only the times in the uncontaminated region
    keepLeft = false;
    keepRight = false;
    if (idxBlock == 1)
      keepLeft = true;
    end
    if (idxBlock == nBlocks)
      keepRight = true;
    end
    eventKeep = true(size(eventTime));
    if ~keepLeft
      eventKeep(eventTime < snipLen) = false;
    end
    if ~keepRight
      eventKeep(eventTime > blockSize-snipLen) = false;
    end
    all_amp{idxBlock} = a(:,eventKeep);
    all_time{idxBlock} = eventTime(eventKeep) + scanNumFrom - 1;
    all_peak{idxBlock} = peakVal(eventKeep);
    all_chan{idxBlock} = peakChan(eventKeep);
    all_pp{idxBlock} = pp(:,:,eventKeep);
    all_err{idxBlock} = err(eventKeep);

    %% Save residual to disk
    if(isSaveRes)
      % Convert to uint16
      waveInt = uint16((wave_res - merecinfo.offset)/merecinfo.slope +...
        merecinfo.minsample);
      if (N < blockSize)
        waveInt = waveInt(:,1:N);
      end
      if ~keepLeft
        waveInt(:,1:snipLen) = [];
      end
      if ~keepRight
        waveInt(:,end-snipLen+1:end) = [];
      end
      count = fwrite(fidRes,waveInt,'uint16');
      if (count < numel(waveInt))
        error('Error writing residual file');
      end
    end
  end

  % Throw out intervals that have no events (to avoid concatenation
  % errors), and consolidate into big matrices
  keepFlag = ~cellfun(@isempty,all_time);
  amplitudes = cat(2,all_amp{keepFlag});
  spiketimes = cat(2,all_time{keepFlag});
  peakVal = cat(2,all_peak{keepFlag});
  peakChan = cat(2,all_chan{keepFlag});
  minmax = cat(3,all_pp{keepFlag});
  err = cat(2,all_err{keepFlag});

  varsToSave = {'fileToFit', 'componentFile', 'optionsForFit', ...
    'amplitudes', 'spiketimes', 'peakVal', 'peakChan', 'minmax', 'err'};
  if isSaveRes
    residualFile = resFile;
    varsToSave(end+1:end+2) = {'residualFile','residual_key'};
  end
  save(fileToSave, '-mat', varsToSave{:});
  if isSaveRes
    fclose(fidRes);
  end
end

function fidRes = ftcmp_openresidual(resFile,header,channels,identifier)
% Update the channel list
  headertext = header.wholeheader;
  headertext = update_value(headertext, 'channel list', num2str(channels));
  % Mark it with an identifier key, so that pairing to fitting parameters
  % can be verified.
  [headertext,succeeded] = update_value(headertext,'identifier key',identifier,false);
  if ~succeeded
    headertext = sprintf('%s\nidentifier key=%d\n',headertext,identifier);
  end
  % Open the residual .merec file and write the header
  [fidRes,msg] = fopen(resFile,'w');
  if (fidRes < 0)
    error(msg);
  end
  count = fwrite(fidRes,headertext,'char');
  if (count < length(headertext))
    error('Couldn''t write the whole header');
  end
  padSize = header.headersize - ftell(fidRes);
  if (padSize > 0)
    count = fwrite(fidRes,zeros(1,padSize),'uint8');
    if (count < padSize)
      error('Couldn''t pad to appropriate size');
    end
  end
end

function options=fillOptions(options)
  options = default(options,'blockSize',2^14);
  %options = default(options,'isUsePeakCriteria',true);
  options = default(options,'isSaveRes',true);
  %options = default(options,'isMergeAdjacent',true);
  %options = default(options,'adjacentNScans',5);
  options = default(options,'fitIterMax',2);
  options = default(options,'debug',false,'maxBlocks',Inf);
  options = default(options,'max_components',Inf);
end
