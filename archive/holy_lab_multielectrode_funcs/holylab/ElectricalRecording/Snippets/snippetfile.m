function [tsnip,snip] = snippetfile(infilename,scanrange,channels,thresh,varargin)
% NOTE - edited 2-13-15 by AM to prevent throwing error on non-Unix
% machines. 
% SNIPPETFILE: snip raw Merec data
% 
% Syntax:  
%    [tsnip,snip] = snippetfile(infilename,scanrange,channels,thresh,varargin)
% 
% pre:
%    infilename: a string (filename), a file identifier (from fopen),
%                or a nchans-by-nscans matrix (which contains the waveform to be
%                processed.) Note snippetfile uses integer arithmetic for
%                speed, so your waveforms should not round to 0! (I.e.,
%                don't convert to volts.)
%    scanrange: a nranges-by-2 matrix, where the waveform between scanrange(i,1) and
%               scanrange(i,2) is analyzed for spikes. Portions of the waveform which do not appear
%               between any of these values are ignored.  Scans are
%               numbered starting at 0.  One may also supply the string
%               'all' to use the entire scanrange.
%    channels: a vector of channel numbers. If waveforms are input, this is
%               an index to the nchans-by-nscans matrix.
%    thresh: a 2-row matrix, each col specifies the boundary values needed to
%            trigger detection of a spike. Row 1 is for troughs and row 2
%            is for peaks.
%    varargin: Variable arguments can be supplied in any order:
%       sniprange: a 2-vector specifying the span of a snippet around the peak. For
%               example, [-10 30] indicates that a snippet starts 10 samples before the 
%               peak and continues until 30 samples after the peak.
%       an options structure: see snipoptions
%       outfilename: a string giving the name of the file to write in saving snippets.
%
% post: (the two return values available only when not sending output to
%        a file)
%    tsnip: rough times for the snippets, in scans. The first scan is 1
%      (i.e., unity-offset).
%    snip: the waveforms of the snippets 
%
% See also: SNIPOPTIONS,SNIPPETFILEAUTO.
  
% History:
%   2004-08-09: 1. Increased buffersize by 2 whenever interpsnips turned on
%               2. Increased minimum bufferize by 4 to allow for rebound
%                  filter
%               (RCH)
%   ?: this file is modified by Jason to handle merec data under Linux
% 
% Note that: 
%       Always use LFS to handle Merec data even the file size < 2G;
%       LFS code works only on Linux;
%       Assume total # of scans < 2G ;
%       Assume snippet file size < 2G ;
%       When processing Merec data, argument infilename must be a file name
%           or fid (if fid, in fact I use fopen() to convert it to file
%           name before call openlfs() )
%           
% limitation: input file must be a Merec data file, and output snippet file
%             will have new ascii format header.
%             If snipping old data is desired, please use snippetfile.m (version
%             200303.200) with WriteSnipHeader, which producing snippet
%             file in old format

  [sniprange,outfilename,soptions] = snippetfileparse(varargin);
  sniprangein = sniprange;
  if (nargout == 0 & ~soptions.tofile)
    warning('No output will be produced, aborting');
    return;
  end
  
  
  %test if the file is merec data
  isMerec=0; % 
  fd=-1; % file descriptor used by LFS in c++ code
  if(isnumeric(infilename) & length(infilename) == 1)
     fseek(infilename,0,'bof');    %back to beginning
     tMagic = fread(infilename, [1,5], 'char');
     if(strcmp(char(tMagic), 'MEREC'))
        isMerec=1; 
        %%% error('Please pass Merec data filename instead of fid to func snippetfile()');
        infilename = fopen(infilename); % convert fid infilename back to file name
     end % if, the file is Merec data
  % end_of_if, fid is passed
  elseif(ischar(infilename))
     if(should_use_lfs(infilename))
        [tFid,tMessage] = openlfs(infilename);
        if(tFid<0) error(tMessage); end
        tMagic = readcharlfs(tFid, 1, [0 4], 0);
        closelfs(tFid);
     else
        [tFid, tMessage]=fopen(infilename);
        if(tFid<0) error(tMessage); end
        tMagic = fread(tFid, [1,5], 'char');
        fclose(tFid);
     end
     if(strcmp(char(tMagic), 'MEREC'))
        isMerec=1; 
     elseif length(soptions.use_alt_header_file) > 1
         isMerec=1;  % Here, the tMagic check is probably failing; let them do whatever recovery they're trying to do w/o inteferance
     end %if, the file is Merec data
  end % elseif, filename is passed
  
  if(isMerec)
     if(~isunix)
%         error('Currently merec data can be processed only on Linux');         %%% !!! commented out by AM 2-15-15 to allow running on Windows 
     end
     fd=openlfs(infilename); % return file descriptor instead of fid
     tInputFileName=infilename;
  end
  
  
  fromfile = 1; % true if read from file (using filename/fid/fd). False if read from mem
  fromfid = 0;  % true if read from file and fid is passed as argument. 
  if (~isMerec & isunix & ~isnumeric(infilename))
    error(['Under Linux, open file first (allow endian conversion) and' ...
           ' pass fid']);
  end
  if (isnumeric(infilename) & length(infilename) == 1)
    fromfid = 1;
    fid = infilename;
    fseek(fid,0,'bof');    % rewind the file
  end
  if (isnumeric(infilename) & length(infilename) > 1)
    fromfile = 0;
    [nchan,nscans] = size(infilename);
    % Make a fake header
    % 
    h.numch = nchan;
    h.nscans = nscans;
    h.channels = 1:nchan;
    wavein = infilename;  % This just creates a reference
  else
      if length(soptions.use_alt_header_file) == 1
        h = readheader(infilename);  
      else
          h = readheader(soptions.use_alt_header_file);
      end
    datastart = h.headersize;
    
    soptions.Fs = h.scanrate;
  end
  
  if (length(channels) ~= size(thresh,2))
    error('Number of channels and number of thresholds do not agree');
  end
  
  % Get the channel index
  [comchan,chanIndex] = intersect(h.channels,channels); % @bookmark: snippetfile_200304081446
  nchannels = length(chanIndex);
  if (nchannels ~= length(channels))
    warning('Not all selected channels were recorded');
  end

  % make items of chanIndex appear increasingly. That is, stable intersect
  % for h.channels
  chanIndex = sort(chanIndex);
  comchan   = h.channels(chanIndex);
  
  tNewThresh= zeros(size(thresh,1), size(comchan,2));
  for i=1:length(comchan) tNewThresh(:,i)=thresh(:,find(channels==comchan(i))); end
  thresh = tNewThresh;
  
  % Set any unsupplied options
  soptions = snipoptions(soptions);
  
  % Check that scanranges are OK
  if ischar(scanrange)
    if strcmp(scanrange,'all')
      scanrange = [0 h.nscans-1];
    else
      error('Unrecognized string for scanrange');
    end
  end
  [nranges,szrange] = size(scanrange);
  if (szrange ~= 2)
    if (nranges == 2)
      scanrange = scanrange';
      nranges = szrange;
    else
      error('scanrange must be a n-by-2 matrix');
    end
  end
  if (max(max(scanrange)) >= h.nscans)
    error('scanrange exceeds size of file');
  end
  
  % Check that snipsize is consistent with other parameters
  snipsize = diff(sniprange) + 1;
  if (snipsize > soptions.blocksize)
    error('blocksize must be bigger than the snip size!');
  end
  
  % Check detection filters, adjust sniprange if necessary
  lendetfilt = 0;
  if ~isempty(soptions.detfilt)
    [lendetfilt,filtchannels] = size(soptions.detfilt);
    %if (lendetfilt ~= snipsize)
    %  error('Detection filter must be same size as snippets');
    %end
    if (filtchannels ~= nchannels)
      error('Number of channels in filter must equal length(channels)');
    end
    % If doing detection filtering, adjust sniprange to compensate for
    % filter
    sniprange = sniprange - lendetfilt;
  end
  
  % If not conditioning, make sure any scaling is applied to thresholds 
  lencondfilt = max(length(soptions.condfilta),length(soptions.condfiltb))-1;
  if (lencondfilt == 0)
    thresh = thresh/(soptions.condfiltb/soptions.condfilta);
  end
  
  % if need transpose:
  if(size(thresh, 1) ~= 2) 
     thresh = thresh';
  end
  
  % To save time & memory, avoid copying & creating new submatrices
  % for data when only a subset of channels are processed---so pass
  % a channel index for working on the raw data, rather than
  % making a submatrix.
  % Because there are several possible sequences of conditioning
  % and detection filtering, we need to set up different channel
  % indices depending on what has already been done
  % wave = indices for "raw" waveform (perhaps after conditioning filtering)
  nchwave = nchannels;
  chIndwave = 1:nchannels;
  if (lencondfilt == 0)
    chIndwave = chanIndex;
    nchwave = h.numch;
  end
  % det = indices for peak detection (perhaps after detection filtering)
  nchdet = nchannels;
  chInddet = 1:nchannels;
  if (lendetfilt == 0)
    chInddet = chIndwave;
    nchdet = nchwave;
  end

  % @jason: chanIndex(1:end) are indices pointing to items in h.channels, 
  % that is, h.channels(chanIndex) is the channel names/IDs (i.e. 0..63)
  % to be cut snippets. @see bookmark snippetfile_200304081446. 
  % so chanIndex(1:end) are also indices pointing to rows of wave.
  % That is, wave(chanIndex(i),:) is the data for channel h.channels(chanIndex(i)).
  % (let j=chanIndex(i), then j-th channel's data is the j-th row)
  %
  % chIndwave and chInddet are similar to chanIndex.
  % 
  % @todo: if chIndwave ~= chInddet, the real channel order of findpeaks()
  % is not same as cutsnippets(). And this will cause problem.
  % 
  % Because findpeaks() and cutsnippets() loop over chInddet and
  % chIndwave, we have to adjust the channel order in snip header:
  tNewChannelList=h.channels(chanIndex);
  if isfield(h,'wholeheader')
    h.labels=split_label(key2value(h.wholeheader, 'label list'));
    tNewLabelList=h.labels(chanIndex, :);
  end
  % then later use tNewLabelList and tNewChannelList to update snip header, 
  % @see bookmark snippetfile_200304081659
  
  % Set up temporary file...
  if (soptions.tofile)
    tname = dirbyname(soptions.tempfilename);
    if ~isempty(tname)
      delete(soptions.tempfilename);
    end
    [fidtemp,message] = fopen(soptions.tempfilename,'w+');
    if (fidtemp == -1)
      warning('May not have permission to write in this directory...');
      error(message);
    end
    nblocks = sum(ceil(diff(scanrange,1,2)/soptions.blocksize));
    tempn = zeros(nblocks,nchannels);
    temppost = zeros(nblocks,nchannels);
    tempposs = temppost;
    tempposFinetime=temppost;
    tempposDetPeak=temppost;
  else
    % ...or set up record in memory
    tsnip = cell(1,nchannels);
    if (nargout > 1)
      snip = cell(1,nchannels);
      for j = 1:nchannels
        snip{j} = zeros(snipsize,0);
      end
    end
  end
  
  % Now the ugly part: have to be able to collect spikes
  % that occur near boundary of a block. (Finite memory limitation.)
  % Buffer the region near boundaries. Set this up now
  
  width = 4;        % sets minimum buffer size (to allow for 3-pt maxima and
                    % rebound spike detection)
  if soptions.interpsnips==1
      bufferleft = max([-sniprange(1)+1,width]);  % if interpsnips is on, you need to make sure you
      bufferright = max([sniprange(2)+1,width]);  % get one more value on either end than you'd otherwise need
  else
      bufferleft = max([-sniprange(1),width]);
      bufferright = max([sniprange(2),width]);
  end
  buffersize = bufferleft+bufferright;
  bufferwave = zeros(nchwave,buffersize);
  bufferdet = zeros(nchdet,buffersize);
  
  % Now get processing! Loop over ranges
  blockindx = 1;
  reboundStats = [];          % variable to collect info on spikes tossed by rebound filter
  for i = 1:nranges
    scanstart = scanrange(i,1);
    if (fromfid)
      fseek(fid,datastart+scanstart*h.numch*2,'bof');
    end
    % Initial conditions for filters
    zicond = zeros(lencondfilt,nchannels);
    zidet = zeros(size(soptions.detfilt) - [1 0]);
    % Process a single contiguous range of time
    while (scanstart + buffersize < scanrange(i,2))
      scanend = min(scanstart + soptions.blocksize - 1,scanrange(i,2));
      % Now read in waveform in block
      % This call can be done natively in new version
      % of matlab, and allows endian conversion
      % (use the fid input style)
      if fromfile
        if ~isunix % not unix:
          wave = readint16(infilename,h.numch,[scanstart scanend]);
        else % is unix:
          if(isMerec)
             wave = readint16lfs(fd,  h.numch,[scanstart scanend], datastart); 
          else
             wave = fread(fid,[h.numch (scanend-scanstart+1)],'*int16');
          end
        end
      else % read from mem:
        wave = int16(wavein(:,scanstart+1:scanend+1));
      end
      % Perform conditioning filtering, if desired 
      if (lencondfilt > 0)
        [wave,zicond] = filterint16(soptions.condfiltb,soptions.condfilta,wave,chanIndex,zicond);
      end
      % Perform detection filtering, if desired
      if isempty(soptions.detfilt)
        wavedet = wave;   % this actually defines a reference, so it's fast
      else
        [wavedet,zidet] = filterint16(soptions.detfilt,[],wave,chIndwave,zidet);
      end
      % Find peaks in the forbidden region from the previous round
      tempdet = [bufferdet,wavedet(:,1:buffersize)];
      [tpeak1,xmax1,detPeaks1,reboundStats] = findpeaks(tempdet,chInddet,thresh,bufferleft, ...
                               bufferright,soptions,reboundStats);
      % Find peaks in the big chunk
      [tpeak2,xmax2,detPeaks2,reboundStats] = findpeaks(wavedet,chInddet,thresh,bufferleft, ...
                               bufferright,soptions,reboundStats);
      % Adjust the scan number indicies of the peaks noted in reboundStats
      % to take into account the scanstart value
      rebounds = [];  
      if size(reboundStats,2) > 0
        [notYetAdjusted] = find(reboundStats(5,:));
        rebounds = reboundStats(:,notYetAdjusted);      % nec. only for testing code below
        reboundStats(4:6,notYetAdjusted) = reboundStats(4:6,notYetAdjusted)+scanstart;
        reboundStats(7,notYetAdjusted) = 0;
      end
      % Cut snippets, if desired
      if (soptions.tofile | nargout > 1)
        tempwave = [bufferwave,wave(:,1:buffersize)];
        sniptemp1 = cutsnippets(tempwave,chIndwave,tpeak1,xmax1,sniprange,soptions);
        sniptemp2 = cutsnippets(wave,chIndwave,tpeak2,xmax2,sniprange,soptions);
      end
      
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%%%%%%% TESTING CODE: To use, comment out last real command in rebound
%       %%%%%%%%%%%%%%%%%%%%%%%% spikes filter section in findpeaks.m, &
%       %%%%%%%%%%%%%%%%%%%%%%%% uncomment this section - then insert stop at the 
%       %%%%%%%%%%%%%%%%%%%%%%%% figure(1) command near the top of the loop and
%       %%%%%%%%%%%%%%%%%%%%%%%% run in debugger mode, using the 'save and run' 
%       %%%%%%%%%%%%%%%%%%%%%%%% button to step through the data one figure's worth
%       %%%%%%%%%%%%%%%%%%%%%%%% at a time.  Can't handle multiple channels.
%       
%       % plot the snippets as they're generated, with an inication of 1) 
%       % where the peak really is in the snippet, 2) whether it
%       % triggered the rebound filter or not, and 3) if it did, the locaiton
%       % of the three points of the triad
%       
%       % step through the plotting 9 at a time
%       
%       toPlot = [sniptemp1{1} sniptemp2{1}];
%       peaks = [tpeak1{1} tpeak2{1}];
%       peaks(2:3,:) = zeros;    %find which peaks triggered rebound filter 
%       for n = 1:size(rebounds,2)
%           triggeredFilter = find(peaks(1,:)==rebounds(6,n));
%           if triggeredFilter > 0
%               peaks(2,triggeredFilter) = 1;
%               peaks(3,triggeredFilter) = n;
%           end
%       end
%       count=0;
%       while count < size(toPlot,2)
%           figure(2)
%           clf
%           for n = 1:9                           % n says which subplot you're on
%               count = count+1;                  % count says which snippet you're on
%               if count <= size(toPlot,2)
%                   subplot(3,3,n)
%                   plot(sniprange(1):sniprange(2),toPlot(:,count),'k')
%                   hold on
%                   yl = ylim;
%                   plot([0 0],[yl(1) yl(2)],'k')
%                   if peaks(2,count) == 1        % if triggered detector
%                       triad = rebounds(4:6,peaks(3,count)) - peaks(1,count)*ones(3,1);
%                       plot([triad(1) triad(1)],[yl(1) yl(2)],'b')
%                       plot([triad(2) triad(2)],[yl(1) yl(2)],'g')
%                       plot([triad(3) triad(3)],[yl(1) yl(2)],'r')
%                       title(['rebound w/ triad points '...
%                               num2str(triad(1)) ', ' num2str(triad(2)) ' and ' num2str(triad(3)) ' (curv=' ...
%                               num2str(rebounds(1,peaks(3,count))) ')'])
%                   elseif size(rebounds,2)>0     % if didn't trigger but may overlap w/some
%                       range = [peaks(1,count)+sniprange(1) peaks(1,count)+sniprange(2)];
%                       withinrange1 = find( (range(1) < rebounds(4,:)) & (rebounds(4,:) < range(2)) );
%                       withinrange2 = find( (range(1) < rebounds(5,:)) & (rebounds(5,:) < range(2)) );
%                       withinrange3 = find( (range(1) < rebounds(6,:)) & (rebounds(6,:) < range(2)) );
%                       for nn = 1:size(withinrange1,2)
%                           x = rebounds(4,withinrange1(nn))-peaks(1,count);
%                           plot([x x],[yl(1) yl(2)],'b')
%                       end
%                       for nn = 1:size(withinrange2,2)
%                           x = rebounds(5,withinrange2(nn))-peaks(1,count);
%                           plot([x x],[yl(1) yl(2)],'g')
%                       end
%                       for nn = 1:size(withinrange3,2)
%                           x = rebounds(6,withinrange3(nn))-peaks(1,count);
%                           plot([x x],[yl(1) yl(2)],'r')
%                       end
%                       title(['peak at: ' num2str(peaks(1,count))])
%                   else                          % in case rebounds empty this round
%                       title(['peak at: ' num2str(peaks(1,count))])
%                   end
%               end
%           end
%       end
%       
%       %%%%%%% END TESTING CODE
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Amalgamate 
      for j = 1:nchannels
          tpeak{j} = [tpeak1{j}-buffersize tpeak2{j}] + scanstart - lendetfilt;
          xmax{j} = [xmax1{j} xmax2{j}];
          detPeaks{j} = [detPeaks1{j} detPeaks2{j}];
          
          % from version 200303.310, don't correct time:
          % if(soptions.interptimes)
          %   for k = 1:length(tpeak{j})
          %   % for k = 1:length(xmax{j})
          %         tpeak{j}(k) = tpeak{j}(k) + xmax{j}(k);
          %   end
          % end
      end
     
      % Append to record, either in memory or in temporary file
      if (soptions.tofile)
        % Do the recordkeeping that will make life easier later,
        % saving file positions & number of snippets on
        % each channel & block
        for j = 1:nchannels
          tempn(blockindx,j) = length(tpeak{j});
          temppost(blockindx,j) = ftell(fidtemp);
          fwrite(fidtemp,tpeak{j},'int32');
          tempposs(blockindx,j) = ftell(fidtemp);
          fwrite(fidtemp,sniptemp1{j},'int16');
          fwrite(fidtemp,sniptemp2{j},'int16');
          
          % more content than old snippet file:
          % +1 (i.e. the 1st additional content): fine times:
          if(soptions.interptimes)
            tempposFinetime(blockindx,j) = ftell(fidtemp);
            fwrite(fidtemp, xmax{j}, 'float32');
          end
          
          % +2: det peaks:
          tempposDetPeak(blockindx,j) = ftell(fidtemp);
          fwrite(fidtemp, detPeaks{j}, 'int16');
          
        end
      else
        % Append to record in memory
        
        for j = 1:nchannels
          tsnip{j}(end+1:end+length(tpeak{j})) = tpeak{j};
          if (nargout > 1)
            snip{j}(:,end+1:end+size(sniptemp1{j},2)) = sniptemp1{j};
            snip{j}(:,end+1:end+size(sniptemp2{j},2)) = sniptemp2{j};
          end
        end
      end
      % Generate the buffer region for next time
      bufferdet = wavedet(:,end-buffersize+1:end);
      bufferwave = wave(:,end-buffersize+1:end);
      scanstart = scanstart + soptions.blocksize;
      blockindx = blockindx + 1;
    end        % current range
  end        % all ranges
  
  % If saving to a file, permute the temporary file
  % to group all snippets on a given channel together
  if soptions.tofile
    disp('Done snippeting, now processing the temporary file...')
    [fidout,message] = fopen(soptions.outfilename,'w');
    if (fidout == -1)
      warning('Can''t open output file. Permission error?');
      error(message);
    end
    % Create and write the header
    snipheader = h;
    snipheader.sniptype = 0;
    snipheader.numofsnips = sum(tempn);
    snipheader.snipbeginoffset = sniprangein(1);
    snipheader.snipendoffset = sniprangein(2);
    snipheader.thresh = thresh;
    snipheader.channels = comchan;
    if ~fromfile
    % if ~isfield(h,'type')
      % The header was created, data was in memory
      warning('Making up header data...');
      snipheader.type = 1;
      snipheader.version = 2;
      %snipheader.channels = -h.channels;
      snipheader.scanrate = 1;
      snipheader.scalemult = 1;
      snipheader.scaleoff = 0;
      snipheader.date = '';
      snipheader.time = '';
      snipheader.usrhdr = '';
    end
    
    tstrHeader=stringize_snip_header(snipheader, soptions, tInputFileName);
    
    % @bookmark: snippetfile_200304081659
    % update header w/ new channel list and label list:
    tstrHeader=update_value(tstrHeader, 'channel list', num2str(tNewChannelList));
    tstrLabel=strcat(tNewLabelList(1,:), ''); % use strcat() to remove
                                              % the tailing zeros
    for i=2:size(tNewLabelList,1)
       % tstrLabel=strcat(tstrLabel, ',', tNewLabelList(i,:));
       tstrLabel=strcat(tstrLabel, '$', tNewLabelList(i,:)); % modified on 20031028
    end
    tstrHeader=update_value(tstrHeader, 'label list', tstrLabel);
    
    
    tTimesfpos=zeros(1, nchannels); % in my notation, t means temp
    tSnipfpos =tTimesfpos;
    if(~soptions.interptimes)
       tFinetimesfpos=[];
    else
       tFinetimesfpos=tTimesfpos;
    end
    tDetpeaksfpos=tTimesfpos;
    
    update_header(fidout, tstrHeader); % write the incomplete header to
                                       % file as place holder
    % fseek(fidout, snipheader.headersize, 'bof');
    %   I had thought use fseek() is enough because in C++ lseek() beyond
    %   end-of-file will append junk to the file and new file position is
    %   correct. But here in matlab, fseek() beyond eof will return
    %   error. So I have to use update_header() to write incomplete header first 
    %   as place holder. Then I don't need fseek() here.
    
    for i = 1:nchannels
      % Write the times
      tTimesfpos(i) = ftell(fidout);
      for j = 1:nblocks
              fseek(fidtemp,temppost(j,i),'bof');
              temptime = fread(fidtemp,tempn(j,i),'int32');
              fwrite(fidout,temptime,'int32'); 
      end
      
      % Write the snippets
      tSnipfpos(i) = ftell(fidout);
      for j = 1:nblocks
              fseek(fidtemp,tempposs(j,i),'bof');
              tempsnip = fread(fidtemp,[snipsize,tempn(j,i)],'int16');
              fwrite(fidout,tempsnip,'int16');
      end
      
      % more content than old snippet file---write fine times and detpeaks:
      % +1 (i.e. the 1st additional content): write fine time:
      if(soptions.interptimes)
         % then write fine time:
         tFinetimesfpos(i) = ftell(fidout);
         for j = 1:nblocks
            fseek(fidtemp,tempposFinetime(j,i),'bof');
            tempFinetime=fread(fidtemp,tempn(j,i),'float32');
            fwrite(fidout,tempFinetime,'float32');
         end
      end
      
      % +2: write det peaks:
      tDetpeaksfpos(i) = ftell(fidout);
      for j = 1:nblocks
         fseek(fidtemp,tempposDetPeak(j,i),'bof');
         tempDetpeak=fread(fidtemp,tempn(j,i),'int16');
         fwrite(fidout,tempDetpeak,'int16');
      end
      
    end % for,
    
    tstrHeader=update_value(tstrHeader, 'timesfpos', num2str(tTimesfpos));
    tstrHeader=update_value(tstrHeader, 'snipsfpos', num2str(tSnipfpos));
    tstrHeader=update_value(tstrHeader, 'finetimesfpos', num2str(tFinetimesfpos));
    tstrHeader=update_value(tstrHeader, 'detpeaksfpos', num2str(tDetpeaksfpos));
    update_header(fidout, tstrHeader); % update the header in file
    
    % Now clean up
    if(isMerec)
       closelfs(fd);    
    end
    fclose(fidtemp);
    fclose(fidout);
    delete(soptions.tempfilename);
    if ischar(infilename)
      disp(['Done with file ',infilename]);
    else
      disp('Done');
    end
  end
  
  % Depending on what options were selected in reboundO, save and/or
  % display any requested feedback
  if (soptions.reboundS ~= 0 & soptions.reboundC ~=0)           % don't bother if screen turned off
      if size(reboundStats,2) < 1
%          print('no rebound spikes found')
          reboundStats = zeros(7,1);
      end
      if (soptions.reboundO == 1 | soptions.reboundO == 3)
          reboundStats = reboundStats(1:4,:);
          save([infilename, '.mat'], 'reboundStats', 'soptions')
      end
      if (soptions.reboundO == 2 | soptions.reboundO == 3)
          [sortt] = find(reboundStats(3,:));
          notPulled = reboundStats(:,sortt);
          pulled = reboundStats;
          pulled(:,sortt) = [];
          nPulled = num2str(size(pulled,2));
          nNotPulled = num2str(size(notPulled,2));
          nCurv = num2str(soptions.reboundC);
          figure
          cmap = [0 0 0;1 0 0];
          colormap(cmap)
          subplot(2,1,1)
          [pulledT,xT] = hist(pulled(2,:),0:(2*soptions.reboundS/100):(2*soptions.reboundS));
          if size(xT,2)<2
              xT = xT';
          end
          [notPulledT] = hist(notPulled(2,:),0:soptions.reboundS/100:soptions.reboundS);
          [zeross] = hist([],0:soptions.reboundS/100:soptions.reboundS);
          h = bar(xT,[pulledT; notPulledT; zeross]','stacked');
          xlim([xT(1) xT(size(xT,2))])
          xlabel('time gap in scan number')
          title([infilename ': Time gap between 1st and 3rd of all triads found (black = ' nPulled ' spike(s) tossed; red = ' nNotPulled ' triad(s) retained)'])
          subplot(2,1,2)
          xMax = 5*soptions.reboundC;
          [pulledC,xC] = hist(pulled(1,:),0:xMax/100:xMax);
          [notPulledC] = hist(notPulled(1,:),0:xMax/100:xMax);
          [zeross] = hist([],0:xMax/100:xMax);
          h = bar(xC,[pulledC; notPulledC; zeross]','stacked');
          xlim([0 xMax])
          xlabel('curvature values (last bar includes all those past axis dimensions, as well)')
          title(['Curvature values of all triads found (black = ' nPulled ' spike(s) tossed; red = ' nNotPulled ' triad(s) retained because curvature greater than threshold = ' nCurv ')'])
      end
  end
  
function [sniprange,outfilename,soptions] = snippetfileparse(P)
  sniprange = [];
  outfilename = '';
  soptions = [];
  inum = zeros(1,length(P));
  istruc = inum;
  ich = inum;
  for i = 1:length(P)
    inum(i) = isnumeric(P{i});
    istruc(i) = isstruct(P{i});
    ich(i) = ischar(P{i});
  end
  srindx = find(inum);if (length(srindx) > 1), error('Can only have one sniprange');end
  if (~isempty(srindx))
    sniprange = P{srindx};
  end
  opindx = find(istruc); if (length(opindx) > 1), error('Can only have one options'); end
  if (~isempty(opindx))
    soptions = P{opindx};
  end
  fnindx = find(ich); if (length(fnindx) > 1), error('Can only have one filename'); end
  if (~isempty(fnindx))
    outfilename = P{fnindx};
    soptions.tofile = 1;
    soptions.outfilename = outfilename;
  end
  return
