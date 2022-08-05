function ssnpFile=fake_channel(fittingResultFile, sortResultDir, options)
% ssnpFile=fake_channel(fittingResultFile)
% PRE:
%    fittingResultFile: a .fit file
%    sortResultDir: the dir to hold the faked autosort result
%    options: a struct has only two fields:
%       onlyGenSsnp=0: this is useful for generate multifile autosort
%                      result. SEE: fake_channel_multifile().
%       nClosestTemplates=5: for each fake channel, #templates to consider
% POST:
%    ssnpFile: the filename of generated .ssnp file

   % fittingResultFile='cycle2_berry.fit'; % for debug purpose
   
   if(nargin==1)
      sortResultDir=''; %
   end
   if(nargin<=2)
      options=[];
   end
   default_options('onlyGenSsnp', 0);
   default_options('nClosestTemplates', 5);
   
   fittingResultVars=load(fittingResultFile, '-mat');
   fileToFit=fittingResultVars.fileToFit; % the .merec file
   templateFile=fittingResultVars.templateFile;
   fitting=fittingResultVars.fitting;
   spikeCollection=fittingResultVars.spikeCollection;
   T=fittingResultVars.T;
   ShiftMatrix=fittingResultVars.ShiftMatrix;
   
   templateVars=load(templateFile, '-mat');
   templates=templateVars.templates;
   channels=templateVars.channels;
   thresh=templateVars.thresh;
   medv=templateVars.medv;

   rawClusterFile=replace_extension(templateFile, '.raw_cluster');
   rawClusterVars=load(rawClusterFile, 'options', '-mat'); % only load variable "options"
   rawClusterOptions=rawClusterVars.options;
  
   nTemplates=length(templates);
   nTemplatesToFake=length(fitting);
   
   % TODO: temp for debug
   % nTemplatesToFake=1;
   
   % gen header
   snipheader=readheader(fileToFit);
   snipheader.sniptype = 0;
   snipheader.numofsnips = zeros(1, nTemplatesToFake);
   snipheader.snipbeginoffset = 1;
   snipheader.snipendoffset = length(channels)*(1+diff(rawClusterOptions.sniprange)); 
   snipheader.thresh = repmat([0 0]', 1, nTemplatesToFake); % TODO: may cause problem
   snipheader.channels = 1:nTemplatesToFake;
   
   soptions.condfilta=0;
   soptions.condfiltb=0;
   soptions.detfilt=[];
   soptions.polarity=0;
   soptions.close=0;
   soptions.troughdepth=0;
   soptions.peaktrough=0;
   soptions.blocksize=0;
   soptions.interptimes=0;
   soptions.interpsnips=0;
   
   tstrHeader=stringize_snip_header(snipheader, soptions, fileToFit);
   
   % change channel list and label list:
   tstrHeader=update_value(tstrHeader, 'channel list', num2str(snipheader.channels));
   tstrLabel=num2str(snipheader.channels, '%d$');
   tstrHeader=update_value(tstrHeader, 'label list', tstrLabel);
   
   % change header size. This is nec. for faked merec file.
   % TODO: hardcoded number temp used here.
   tstrHeader=update_value(tstrHeader, 'header size', num2str(4*1024*1024));
   
   ssnpFile=replace_extension(fittingResultFile, '_fake.ssnp');
   [fidout,message] = fopen(ssnpFile,'w');
   if (fidout == -1)
      warning('Can''t open output file. Permission error?');
      error(message);
   end
   update_header(fidout, tstrHeader); % write the incomplete header to file as place holder
   
   fakeChannelInfoFile=replace_extension(ssnpFile, '.fake_info');

   if(isempty(sortResultDir))
      sortResultDir='sort_fake'; % TODO:
   end
   
   if(~options.onlyGenSsnp)
      mkdir(sortResultDir); 
   end

   % calc template peak positions
   templatePeakPos=nan(1, nTemplates);
   for templateID=1:nTemplates
      [tt, templatePeakPos(templateID)]=max(abs(templates{templateID}));
   end
   
   
   memm = merecmm(fileToFit,'tovolts',false,'contiguous',true);
   sniprange=rawClusterOptions.sniprange;

   fakeInfo=struct('templateID', {}, 'closestTemplateIDs', {}, ...
      'snipLabels', {} );
   
   for templateID=1:nTemplatesToFake
      fprintf('working on template %d of %d\n', templateID, nTemplatesToFake);
   
      templatesMat=cell2mat(templates);
      curTemplate=templates{templateID};
      curTemplate=repmat(curTemplate, 1, nTemplates);
      distance=sum((templatesMat-curTemplate).^2).^0.5;
      [tt, closestTemplateIDs]=sort(distance);
      if length(closestTemplateIDs)>options.nClosestTemplates
        closestTemplateIDs=closestTemplateIDs(1:options.nClosestTemplates); % only pick closest 4.
      end
      % NOTE: closestTemplateIDs are the template IDs (indices)
      % NOTE: the first is cur template itself
      
      peakPositions=cell(1, length(closestTemplateIDs));
      snipLabels   =cell(1, length(closestTemplateIDs));
      for idx=1:length(peakPositions)
         % eventTimes=       [fitting{closestTemplateIDs(idx)}.eventTime];
         % shifts    =double([fitting{closestTemplateIDs(idx)}.shift] );
         % peakPositions{idx}=eventTimes+shifts+sniprange(1); 
         peakPositions{idx}=[fitting{closestTemplateIDs(idx)}.shiftedTime];
         
         snipLabels{idx}=ones(1,length(peakPositions{idx}))*idx; 
      end % for, each template to compare

      % NOTE: the time is not sorted. 
      % TODO: maybe we should sort the time and re-arrange the snippets
      peakPositions=cat(2, peakPositions{:});
      snipLabels=cat(2, snipLabels{:}); % NOTE: snipLabels are the indices to closestTemplateIDs
      
      % NOTE: closestTemplateIDs(snipLabels(i)) is the template
      %     ID for the i-th snippet.
      %       peakPositions(i) is the peak position of the i-th snippet
      snipLen=(diff(sniprange)+1)*length(channels);
      nSnips=length(peakPositions);
      snips=zeros(snipLen, nSnips, 'int16'); % TODO: hardcoded 'int16'
      peakHeights=nan(1, nSnips);
      for idxSnip=1:nSnips
         ttSnip=memm(channels, peakPositions(idxSnip)+sniprange);
         ttSnip=ttSnip';
         snips(:,idxSnip)=ttSnip(:);
         
         % calc detpeak
         peakHeights(idxSnip)=ttSnip(templatePeakPos( closestTemplateIDs(snipLabels(idxSnip)) ));
      end % for, each snippet

      % the times
      finetimes=zeros(1, nSnips, 'single');
      
      % finally we can write the fake data now 
      tstrHeader=snipfile_append_channel(fidout, tstrHeader, templateID, peakPositions, ...
         snips, finetimes, peakHeights);
      
      % save the fake info in mem:
      fakeInfo(end+1)=struct('templateID', templateID, 'closestTemplateIDs', closestTemplateIDs, ...
         'snipLabels', snipLabels );
   
      if(options.onlyGenSsnp)
         continue; % this skips the generating sort_info part
      end
      
      % create channel specific dir
      chanSpecificDir=fullfile(sortResultDir, ['chan' num2str(templateID)]);
      mkdir(chanSpecificDir);
      
      % prepare the data in autosort_info.mat
      sort_info.channel=templateID;
      sort_info.sniprange=[snipheader.snipbeginoffset snipheader.snipendoffset];
      sort_info.use_projection=1;
      sort_info.t2V=1;
      sort_info.landmarkWaveform=[templates{closestTemplateIDs}];
      sort_info.landmarkT=zeros(1, length(closestTemplateIDs));
      sort_info.Rfactor=nan;
      sort_info.landmarkClust=1:length(closestTemplateIDs); % NOTE: this is the clustering assignment to each landmark
      sort_info.timeMarker=[]; % TODO: as a placeholder
      % calc projection directions
      % use proj_with_svd() instead of lda(). Also use
      % nClostestTemplates of prj dir instead of 2.
      
      %[eigenvec, eigenval]=lda(sort_info.landmarkWaveform, sort_info.landmarkClust);
      % sort_info.projectDirections=eigenvec(:, 1:2);
      [eigenvec]=proj_with_svd(sort_info.landmarkWaveform);
      sort_info.projectDirections=eigenvec;
      
      % save sort_info
      sortInfoFilename=fullfile(chanSpecificDir, 'autosort_info.mat');
      save(sortInfoFilename, 'sort_info', '-mat');
      
   end % for, each fake channel
   
   fclose(fidout);

   % save the fake info
   save(fakeChannelInfoFile, 'fakeInfo', 'templateFile', '-mat');

   if(options.onlyGenSsnp)
      return; % this skips the generating overview part
   end
   
   % now save the overview
   overview.sorthead = snipfile2sortheader(ssnpFile);
   overview.options = []; % TODO: as a place holder
   overview.isFake  = 1;
   
   overviewFilename=fullfile(sortResultDir, 'overview.mat');
   save(overviewFilename, '-struct', 'overview', '-mat');
   
   
   
