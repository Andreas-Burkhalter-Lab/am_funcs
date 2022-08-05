function verify_fit(fittingResultFile, options)
% verify_fit(fittingResultFile)
% PRE:
%    fittingResultFile: a *.fit file
   if(nargin==1)
      options=struct;
   end
   options=fillOptions(options);


   fittingResultVars=load(fittingResultFile, '-mat');
   fileToFit=fittingResultVars.fileToFit;
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
   
   
   memm = merecmm(fileToFit,'tovolts',true,'contiguous',true);
   sniprange=rawClusterOptions.sniprange;
   snipMedian=repmat(medv, 1, diff(sniprange)+1);
   %wave=memm(channels, [scanNumFrom scanNumTo]);
   
   header=memm.header;
   scanRangeToVerify=options.timeRangeToVerify*header.scanrate; 
   scanRangeToVerify(1)=scanRangeToVerify(1)+1; % +1: matlab is 1-base numbering.
   if(scanRangeToVerify(1)>header.nscans)
      error('tried to verify fitting beyond file end');
   end
   if(scanRangeToVerify(2)<0 || scanRangeToVerify(2)>header.nscans) 
      scanRangeToVerify(2)=header.nscans;
   end
   
   allEventTimes=unique([spikeCollection.eventTime]);
   eventsIndicesToSee=find(allEventTimes>scanRangeToVerify(1) & allEventTimes<scanRangeToVerify(2));
   
   if(options.type==1)
      figure;
      nSaw=0;
      for eventIdx=eventsIndicesToSee
         if(eventIdx<options.startIndex)
            continue;
         end
         
         % clear snip;
         nSaw=nSaw+1;
         spikes=spikeCollection(eventIdx);
         % shiftedTemplates=T(spikes.rowIdxToShiftMatrix,:)';
         clear shiftedTemplates;
         for tIdx=1:length(spikes.rowIdxToShiftMatrix)
            shiftedTemplates{tIdx}=T(:,:, spikes.rowIdxToShiftMatrix(tIdx))';
            shiftedTemplates{tIdx}=shiftedTemplates{tIdx}(:);
         end
         shiftedTemplates=cat(2, shiftedTemplates{:});
         %shiftedTemplates=T(:,:, spikes.rowIdxToShiftMatrix)';
         %shiftedTemplates=shiftedTemplates(:);
         eventPos=spikes.eventTime(1); % TODO: may off-by-one?
         
         %if(spikes.eventTime(1)<scanRangeToVerify(1) || spikes.eventTime(1)>scanRangeToVerify(2))
         %   continue;
         %end % if, out of the range interested
         
         ttSnip=memm(channels, eventPos+sniprange);
         ttSnip=ttSnip-snipMedian;  % normalize
         ttSnip=ttSnip';
         % snip(:,1)=ttSnip(:);% each column is a snippet; each snippet is saved channel-by-channel
         hLines=plot([ttSnip(:) shiftedTemplates sum(shiftedTemplates,2)]);
         set(hLines(1), 'color', 'blue'); % real data is blue
         set(hLines(end), 'color', 'black', 'lineWidth', 2);  % template is black
         tTitle=sprintf('@%d, event %d of %d: %d template(s) --- [%s], shifts --- [%s]', ...
            spikes.eventTime(1), nSaw, ...
            length(eventsIndicesToSee), ...
            size(shiftedTemplates, 2), ...
            num2str(spikes.templateIdx), num2str(spikes.shift, '%+4d'));
         title(tTitle);
         
         pause;
      end % for, each event
   else
      nTemplates=length(fitting);
      for templateIdx=1:nTemplates
         template=templates{templateIdx};
         clear snip;
         figure;
         
         allSpikeTimes=[fitting{templateIdx}.eventTime];
         spikeIndicesToSee=find(allSpikeTimes>scanRangeToVerify(1) & allSpikeTimes<scanRangeToVerify(2));
         nSaw=0;
         for idxSpike=spikeIndicesToSee
            if(idxSpike<options.startIndex)
               continue;
            end
            
            nSaw=nSaw+1;
            curSpike=fitting{templateIdx}(idxSpike);
            % peakPos=curSpike.shiftedTime-sniprange(1)+1; % NOTE: shiftedTime is left edge's position
            peakPos=curSpike.shiftedTime; % TODO: may off-by-one?
            % peakPos=curSpike.shiftedTime;
            
            %if(curSpike.shiftedTime<scanRangeToVerify(1) || curSpike.shiftedTime>scanRangeToVerify(2))
            %   if(curSpike.shiftedTime>scanRangeToVerify(2))
            %      break;
            %   end
            %   continue;
            %end % if, out of range interested
            
            ttSnip=memm(channels, peakPos+sniprange);
            ttSnip=ttSnip-snipMedian;  % normalize
            ttSnip=ttSnip';
            snip(:,idxSpike)=ttSnip(:);% each column is a snippet; each snippet is saved channel-by-channel
            hLines=plot([template  snip(:,idxSpike)]);
            set(hLines(1), 'color', 'black', 'linewidth', 2); % template is black
            set(hLines(2), 'color', 'blue');  % real data is blue
            
            tTitle=sprintf('@%d, template %d%s, event %d of %d', ...
               curSpike.eventTime, templateIdx, num2str(curSpike.shift, '%+4d'), ...
               nSaw, length(spikeIndicesToSee) ...
               );
            title(tTitle);
            
            pause
         end
         
         % TODO: plot everything on one graph
      end % for, each template
   end % else, 

function options=fillOptions(options)
   default_options('type', 1);
   default_options('timeRangeToVerify', [0 -1]);
   default_options('startIndex', -1);
   
   