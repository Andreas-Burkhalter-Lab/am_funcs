function cluster_events_fine(filename, options)

% cluster_events_fine(filename, options)
%
% PRE: 
%    filename: a *.raw_cluster file or *.fine_cluster file (for resuming).
%
%    options: (optional) a struct that has fields:
%       whatToPca='aboveThresh': use above threshold section to do pca,
%                                otherwise use specified channels to do pca.
%       window_suppress_on: a 0 (default) or a 1, indicating whether you
%                                (a) the window allowing additional
%                                alignment peaks to never be shown, 
%                                and (b) the window that comes up after the
%                                interactive clustering window
%                                allowing the deletion of individual
%                                snippets on the basis of their full
%                                channel spectrum to only be shown in cases
%                                where there are relatively few snippets 
%       abs_pathname_correct_on: a 0 (default) or 1, indicating whether you
%                                want autodetection and correction of
%                                situations where the relative pathnames of
%                                your directories have been constant, but
%                                the absolute pathnames have changed (eg,
%                                the entire experiment has been coppied en
%                                masse to a new sata drive)
%       win_locations:  a structure array of 'Position' specifications for
%                                windows that pop up during clustering; for
%                                an example of their use, please see rclust
%                                (a function that runs cluster_events_fine
%                                with specified options):
%                                   - win_locations.snipToPca 
%                                   - win_locations.PC
%                                   - win_locations.snip
%                                   - win_locations.mdexplore
% HISTORY:
%   ?           (ZG)    wrote it
%   2007-03-27  (RCH)   add option to suppress some of the lesser used
%                       window popups
%   2007-04-13  (RCH)   added abs_pathname_correct_on and win_locations as
%                       options
%   2007-04-20  (RCH)   changed range for snippet display 
%                       (to change again, search for " ranges(:,1) =
%                       overx(ranges(:,1),1); "

   if(nargin==1)
      options=struct;
   end
   options=fillOptions(options);
   fineClusterOptions=options;
   
   [tt,tt,ext] = fileparts(filename);
   if(strncmp(ext, '.raw_cluster',12))
     % We're starting fresh
      rawClusterFilename=filename;
      fileToSave=replace_extension(rawClusterFilename, '.fine_cluster');

      % these 3 fields are changed cross iterations
      rawClusterStartIndex=1;
      templates={};
      fineClusters=struct('idxRawCluster', {}, 'idxFineCluster', {}, ...
         'channelsUsedToAlign', {}, 'alignmentShifts', {}, ...
         'snipIndices', {});
   else
     % We're resuming a previous operation
      fileToSave=filename;
      tt=load(fileToSave, '-mat');
      rawClusterFilename=tt.rawClusterFilename;
      
      rawClusterStartIndex=tt.rawClusterStartIndex;
      templates=tt.templates;
      fineClusters=tt.fineClusters;
   end

   % create the dir where to save snippet files
   snippetDir=replace_extension(fileToSave, '.snippets');
   [tSuc, tMsg, tMsgId] = mkdir(snippetDir);
   if(~tSuc)
      error(['failed to create dir ' snippetDir]);
   end
   
   tt=load(rawClusterFilename, '-mat');
   % snippets=tt.snippets;
   channels=tt.channels;
   if options.abs_pathname_correct_on
       merecFilename = [pwd filesep rawClusterFilename(1:(end-18)) '.merec'];
   else
       merecFilename=tt.merecFilename;
   end
   rawClusterOptions=tt.options;
   labels=tt.labels;
   berryLabels=tt.berryLabels;
   if isfield(tt,'sniptimes')
     sniptimes=tt.sniptimes;
   else
     sniptimes = tt.eventTime;
   end
   thresh=tt.thresh;
   medv=tt.medv;

   have_landmarks = false;
   if isfield(tt,'landmarks')
     landmarks=tt.landmarks;
     landmarkClusters=tt.landmarkClusters;
     landmarkPos=tt.landmarkPos;
     landmarkFinalPos=tt.landmarkFinalPos; %
     have_landmarks = true;
   else
     eventMean = tt.eventMean;
   end
   
   extractSnipOptions=rawClusterOptions;
   extractSnipOptions.maxShift=3; % TODO: make it an option
   
   nRawClusters=length(sniptimes);
   
   if(rawClusterStartIndex>nRawClusters)
      fprintf('fine-clustering has already been done.\n');
      return;
   end
   
   ttStartIndex=rawClusterStartIndex;
   for idxRawCluster=ttStartIndex:nRawClusters
   %for idxRawCluster=3:3    
      fprintf('Loading %5d snippets in raw cluster %3d, please wait ...', ...
         length(sniptimes{idxRawCluster}), idxRawCluster);
      
      %snip=snippets{idxRawCluster};
      snip=extract_snip(merecFilename, channels, idxRawCluster, sniptimes, medv, extractSnipOptions);
      % spikeTimes=sniptimes{idxRawCluster};
      
      fprintf(' done\n');
      
      % skip raw cluster that has too few snippets
      minNumSnippetsInRawCluster=5; % TODO: make it an options
      nSnippets=size(snip,2);
      % fprintf('#snippets in raw cluster %d: %d\n', idxRawCluster, nSnippets);
      if(nSnippets<minNumSnippetsInRawCluster)
         warning(['skip raw cluster ' num2str(idxRawCluster) ' b/c it has only ' ...
            num2str(nSnippets) ' snippets']);
         continue;
      end

      if have_landmarks
        if options.window_suppress_on
          fprintf('...on raw cluster %d of %d; op.window_suppress_on is turned on...',idxRawCluster,nRawClusters)
        else
          figOverlayedSnip=figure;
          plot(snip); % NOTE: with the extra samples for alignment
          % plot landmarks' initial positions
          figLandmarkPos=figure;
          figureTitle=sprintf('Please pick channels for raw cluster %d of %d', idxRawCluster, nRawClusters);
          set(figLandmarkPos, 'name', figureTitle);
          plot(landmarkPos(:, landmarkClusters{idxRawCluster}));

          % plot landmarks' final positions
          hold on;
          plot(landmarkFinalPos(:, landmarkClusters{idxRawCluster}), 'k', ...
            'linewidth', 3); % use black color for final positions
        end

        % calc preset channels for PCA
        averageFinalPos=mean(landmarkFinalPos(:, landmarkClusters{idxRawCluster}), 2);
        presetChannels=calcPresetChannels(averageFinalPos); % NOTE: presetChannels is channel indices, not real channel nums
        

        % if(~isequal(options.whatToPca, 'aboveThresh'))
        if ~options.window_suppress_on
          channelsToSubcluster=pickChannels(figLandmarkPos, presetChannels);
          channelsToSubcluster=unique(channelsToSubcluster);
          delete(figLandmarkPos);
        else
          channelsToSubcluster=presetChannels;
        end

        channelsUsedToAlign=channelsToSubcluster; %
      else
        [tmp,channelsUsedToAlign] = max(abs(eventMean(:,idxRawCluster)));
        channelsToSubcluster = find(abs(eventMean(:,idxRawCluster)) > thresh);
        if isempty(channelsToSubcluster)
          channelsToSubcluster = channelsUsedToAlign;
        end
      end
        
      % align on picked channels only
      snipLen=diff(extractSnipOptions.sniprange)+1+extractSnipOptions.maxShift*2;
      rows=[];
      for ch=make_vector(channelsUsedToAlign, 'row')
         rows=[rows ch*snipLen-snipLen+1:ch*snipLen];
      end
      if ~options.window_suppress_on
        figSnipToAlign=figure;
        set(figSnipToAlign, 'name', 'Unaligned snippets');
        plot(snip(rows, :));
        %startPositions=snipLen:snipLen:(length(rows)-1);
        %plotSeparationLines(gca, startPositions);
        peakPositions = -extractSnipOptions.sniprange(1)+1+extractSnipOptions.maxShift:snipLen: ...
          (length(rows)-1);
        plotSeparationLines(gca,peakPositions);
      end
      shifts=align_snip(snip(rows,:), snipLen, extractSnipOptions.maxShift);
      
      % shift & resample snippets of all channels according to shifts
      snip=resample_snip(snip, shifts, extractSnipOptions);
      
      % now snip is shifted/resampled, resized (each snip is of size
      %    extractSnipOptions.sniprange)
      
      snipLen=diff(extractSnipOptions.sniprange)+1;
         
      % if(~isequal(options.whatToPca, 'aboveThresh'))   
      rows=[];
      for ch=make_vector(channelsToSubcluster, 'row')
         rows=[rows ch*snipLen-snipLen+1:ch*snipLen];
      end
      
      startPositions=snipLen:snipLen:(length(rows)-1);
      
      if(isequal(options.whatToPca, 'aboveThresh'))
         % if use data above threshold to do pca
         
         
         % remove the extra
         % TODO:
         % shifts=zeros(1, size(snip,2));
         % snip=resample_snip(snip, shifts, extractSnipOptions);
         meanSnippet=mean(snip,2);
         tThresh=repmat(thresh, 1, snipLen); 
         tThresh=tThresh';
         tThresh=tThresh(:);
         rows=find(abs(meanSnippet)>=tThresh);
         
         if(length(rows)<3)
            % pick the top 3
            absMeanSnippet=sort(abs(meanSnippet));
            big3=absMeanSnippet(end-2:end);
            peaks=findainb(big3, abs(meanSnippet));
            peaks=sort(peaks);
            rows=[];
            startPositions=[1];
            for peakpos=make_vector(peaks, 'row')
               newRows=peakpos+[-10:18]; % the range is peakpos+[-10 18]
               % remove nonpositive indices here
               newRows(newRows<1)=[];
               nOverlaps=length(intersect(rows, newRows));
               newRows=newRows(nOverlaps+1:end);
               if(isempty(newRows))
                  continue;
               end
               
               startPositions(end)=startPositions(end)-nOverlaps;
               rows=[rows newRows];
               startPositions(end+1)=length(rows)+1;
            end % for each peak
         else
            % rows=sorts(rows); % not necessary b/c find()
            ranges=segment_indices(rows);
            ranges(:,2)=ranges(:,2)+ 18; % 18 extra for under-thresh part
            ranges(:,1)=ranges(:,1)- 10; % add 10 points before the peak, too
            ranges(:,1) = overx(ranges(:,1),1); % make sure hasn't pushed into negative indicies
            rows=[];
            startPositions=[1];
            for range=ranges'
               newRows=range(1):range(2);
               newRows(newRows>length(meanSnippet))=[]; % remove out of range indices
               nOverlaps=length(intersect(rows, newRows));
               newRows=newRows(nOverlaps+1:end);
               if(isempty(newRows))
                  continue;
               end
               
               startPositions(end)=startPositions(end)-nOverlaps;
               rows=[rows newRows];
               startPositions(end+1)=length(rows)+1;
            end % for, each range
         end
         startPositions([1,end])=[]; % remove the 1st and the beyond-the-end
      end % if, use above_thresh method

      
      % plot snippets used to do PCA
      figSnipToPca=figure; 
      if isstruct(options.win_locations)
          set(figSnipToPca,'Position',options.win_locations.snipToPca)
      end
      plot(snip(rows, :));
      set(figSnipToPca, 'name', 'Snippets used in PCA');
      plotSeparationLines(gca, startPositions);
      
      
      [subclusters, aborted]=subcluster_events(snip(rows, :), startPositions, options);
      
      % close figures
      if ~options.window_suppress_on && have_landmarks
          free(figOverlayedSnip)
      end
      free(figSnipToPca);
      if(isequal(options.whatToPca, 'aboveThresh')) & ~options.window_suppress_on
         free(figSnipToAlign);
      end
            
      if(aborted)
         return;
      end
      
      for idxSubCluster=1:length(subclusters)
         curFineCluster.idxRawCluster=idxRawCluster;
         curFineCluster.idxFineCluster=idxSubCluster;
         curFineCluster.channelsUsedToAlign=channelsUsedToAlign;
         curFineCluster.alignmentShifts=shifts; % NOTE: the shifts are for all snippets, not just snippets that belong to this fine cluster
         curFineCluster.snipIndices=subclusters{idxSubCluster}; % NOTE: the indices is valid to the raw cluster only
         curFineClusterSnippets=snip(:, curFineCluster.snipIndices);
         
         % user may want to check data on all channels and get rid of some snippets
         % NOTE: the return value is column indices to curFineClusterSnippets
         if ~options.useMedianToClean
           if (size(curFineClusterSnippets,2)<15) | ~options.window_suppress_on
             [indicesLeft, aborted]=cleanFineCluster(curFineClusterSnippets);
             if(aborted)
               return;
             end
             curFineCluster.snipIndices=curFineCluster.snipIndices(indicesLeft);
             curFineClusterSnippets=curFineClusterSnippets(:,indicesLeft);
           end
         end
         
         % save the snippets
         curSnipFilename=fullfile(snippetDir, genSnipFilename(idxRawCluster, idxSubCluster));
         save(curSnipFilename, 'curFineClusterSnippets', ...
            ... % 'curFineCluster', ...  % duplicate info may cause confusion
            '-mat');
         
         % generate template by averaging the snippets (each fine cluster
         %   generates one template):
         if options.useMedianToClean
           templates{end+1} = median(curFineClusterSnippets,2);
         else
           templates{end+1}=mean( curFineClusterSnippets ,2);
         end
         
         fineClusters(end+1)=curFineCluster;
      end % for, each newly got fine cluster
      
      rawClusterStartIndex=idxRawCluster+1;
      
      % save it to file
      sniprange = extractSnipOptions.sniprange;
      save(fileToSave, 'templates', ...
         'fineClusters', ...
         'rawClusterStartIndex', ...
         'channels', ...
         'rawClusterFilename', ...
         'thresh', ...
         'medv', ...
         'fineClusterOptions', ...
	 'sniprange', ...
         '-mat');
      
   end % for, each raw cluster

   
function [indicesLeft, aborted]=cleanFineCluster(snippets)
   figSnip=figure;
   hLines=plot(snippets);
   axesSnip=gca;
   hold(axesSnip, 'on');
   hAvgLine=plot(mean(snippets,2));
   set(hAvgLine, 'linewidth', 2, 'color', 'black', 'hittest', 'off');
   % for idxLine=1:length(hLines)
   %    set(hLines(idxLine), 'tag', idxLine);
   % end % for, each line
   
   title(['all channels'' data for current template']);
   
   set(hLines, 'ButtonDownFcn', @onClickOnSnip);
   set(axesSnip, 'ButtonDownFcn', @onSelectMultipleSnip);
   
   undoStack={};
   
   % set appdata then waitfor ctrl+q
   setappdata(figSnip, 'undoStack', undoStack);
   setappdata(figSnip, 'hLines', hLines);
   setappdata(figSnip, 'hAvgLine', hAvgLine);
   setappdata(figSnip, 'snippets', snippets);
   
   bind_shortcut(figSnip, 'ctrl+q', @onDone);
   bind_shortcut(figSnip, 'ctrl+a', @onDone);
   bind_shortcut(figSnip, 'delete', @onRemoveSnips);
   bind_shortcut(figSnip, 'ctrl+d', @onRemoveSnips);
   bind_shortcut(figSnip, 'ctrl+z', @onUndoRemoveSnips);
   bind_shortcut(figSnip, 'ctrl+u', @onUnselectAllSnips);
   helpMsg=sprintf(['Press ctrl+q or ctrl+a to finish;\n' ...
      'Press DEL or ctrl+d to remove snippets from the fine cluster;\n' ...
      'Press ctrl+z to undo the deletion;\n' ...
      'Press ctrl+u to unselect all snippets']);
   bind_shortcut(figSnip, 'f1', {@show_help, 'fine tune the template', helpMsg});
   
   set(figSnip, 'UserData', 0);
   waitfor(figSnip,'UserData');
   
   if(ishandle(figSnip))
      % if user press ctrl+q
      indicesLeft=find(isvisible(hLines));
      aborted=false;
   else
      % else user click "X" on title bar, which means user wants to stop
      indicesLeft=[];
      aborted=true;
   end
   
   free(figSnip);

   
function toggleSelection(hLine)
   if(strcmp(get(hLine, 'marker'), '*'))
      set(hLine, 'marker', 'none');
   else
      set(hLine, 'marker', '*');
   end
      
      
function onClickOnSnip(sender, event)
   hLine=sender;
   toggleSelection(hLine);

   
function onSelectMultipleSnip(sender, event)
   axesSnip=sender;
   figSnip=get_parent_fig(axesSnip);
   rect = GetSelRect;
   % loop over all lines and test if a line is in the rect
   hLines=getappdata(figSnip, 'hLines'); % OR: findobj() can be used.
   hLines=hLines(ishandle(hLines));
   for idxLine=1:length(hLines)
      hLine=hLines(idxLine);
      if(is_line_in_rect(hLine, rect))
         toggleSelection(hLine);
      end
   end % for, each line


 
function onUnselectAllSnips(sender, event)
   figSnip=sender;
   hLines=getappdata(figSnip, 'hLines'); % all line objs
   set(hLines, 'marker', 'none');
   
   
   
function onUndoRemoveSnips(sender, event)
   figSnip=sender;
   hLines=getappdata(figSnip, 'hLines'); % all line objs
   
   undoStack=getappdata(figSnip, 'undoStack');
   if(isempty(undoStack))
      uiwait(msgbox('No undo is available','cluster_event_fine','modal'));
      return;
   end
   
   set(undoStack{end}, 'visible', 'on');
   undoStack(end)=[];
   setappdata(figSnip, 'undoStack', undoStack);
   
   % update the avg line
   hAvgLine=getappdata(figSnip, 'hAvgLine');
   snippets=getappdata(figSnip, 'snippets');
   bVisible=isvisible(hLines);
   avg=mean(snippets(:, bVisible), 2);
   set(hAvgLine, 'ydata', avg);
   if(~any(bVisible))
      uiwait(msgbox('All snippets are removed, you may consider undo','cluster_event_fine','modal'));
   end
     
   
function onRemoveSnips(sender, event)
   figSnip=sender;
   hLines=getappdata(figSnip, 'hLines'); % all line objs
   % hLines=hLines(ishandle(hLines));
   hSelectedLines=findobj(figSnip, 'type', 'line', 'Marker', '*', 'visible', 'on'); % find selected lines
   % free(hSelectedLines); 
   set(hSelectedLines, 'visible', 'off');
   
   undoStack=getappdata(figSnip, 'undoStack');
   undoStack{end+1}=hSelectedLines;
   setappdata(figSnip, 'undoStack', undoStack);
   
   % update the avg line
   hAvgLine=getappdata(figSnip, 'hAvgLine');
   snippets=getappdata(figSnip, 'snippets');
   bVisible=isvisible(hLines);
   avg=mean(snippets(:, bVisible), 2);
   set(hAvgLine, 'ydata', avg);
   if(~any(bVisible))
      uiwait(msgbox('All snippets are removed, you may consider undo','cluster_event_fine','modal'));
   end
   
   
  
function filename=genSnipFilename(idxRawCluster, idxSubCluster)
    filename=sprintf('r%04d_f%04d.mat', idxRawCluster, idxSubCluster);
   
% TODO: move this func out of this file
% 
function shifts=align_snip(snip, snipLen, maxShift)
% snip: each col is a snippet; a snippet is saved channel-by-channel; 
%       each channel's snippet has length snipLen
% maxShift: is also the extra # samples for alignment
% shifts: the shift for each column (i.e. for each snippet)
   nSnip=size(snip,2);
   snippets=cell(1, nSnip);
   for idx=1:nSnip
      snippets{idx}=reshape(snip(:,idx), snipLen, []); % in snippets{i}, each col is for each channel
   end
   shifts=zeros(1,nSnip);
   coreSnipLen=snipLen-maxShift*2;
   nChannels=size(snippets{1},2);
   for ite=1:maxShift
      % calc the mean snippet
      tsum=zeros(coreSnipLen, nChannels);
      for idxSnip=1:nSnip
         tsum=tsum+snippets{idxSnip}(shifts(idxSnip)+maxShift+[1:coreSnipLen], :);
      end % for, each snippet
      tmean=tsum/nSnip;
      
      bestRelativeShifts=zeros(1,nSnip);
      for idxSnip=1:nSnip
         errors=zeros(1,3); % 3 is length(-1:1)
         for relativeShift=-1:1
            tSnip =snippets{idxSnip}(shifts(idxSnip)+maxShift+[1:coreSnipLen]+relativeShift, :);
            errors(relativeShift+2)=sum(sum((tSnip-tmean).^2));
         end % for, each relative shift from -1/0/+1
         % TODO: may need fractional shifts (by 3 pt interplation of -1/0/1)
         [tt, idxSmallest]=min(errors);
         bestRelativeShifts(idxSnip)=idxSmallest(1)-2;
      end % for, each snippet
      shifts=shifts+bestRelativeShifts;
      if(~any(bestRelativeShifts))
         break;
      end % if, no shift change
   end % for, each iteration

   
% shift & resample & resize snippets   
% TODO: when use fractional shifts, this func need big change
function snip=resample_snip(origSnip, shifts, extractSnipOptions)
   maxShift=extractSnipOptions.maxShift;
   nSnippets=size(origSnip,2);
   coreSnipLen=diff(extractSnipOptions.sniprange)+1;
   origSnipLen=coreSnipLen+maxShift*2;
   nChannels=size(origSnip,1)/origSnipLen;
   snip=zeros(coreSnipLen*nChannels, nSnippets);
   for idxSnip=1:nSnippets
      origSnippet=reshape(origSnip(:,idxSnip), origSnipLen, []);
      curSnippet=origSnippet(shifts(idxSnip)+maxShift+[1:coreSnipLen], :);
      snip(:,idxSnip)=curSnippet(:);
   end % for, each snippet


function presetChannels=calcPresetChannels(averageFinalPos)
   isUseTop1Only=1; % TODO: make it an option

   if(isUseTop1Only)
      [tmax, tidx]=sort(abs(averageFinalPos));
      presetChannels=tidx(end-0:end); % if top 2: end-1:end 
      presetChannels=sort(presetChannels);
   else
      s=std(averageFinalPos);
      presetChannels=find(abs(averageFinalPos)>s); % NOTE: averageFinalPos is 0-aligned
      % find the peaks
      t=abs(averageFinalPos);
      isPeak=t([1 1:end-1])<t & t>t([2:end end]);
      presetChannels=intersect(presetChannels, find(isPeak));
   end


function channelsToSubcluster=pickChannels(fig, presetChannels)
   hax=gca; % TODO: should use a more robust way
   set(hax, 'ButtonDownFcn', @onAxesBtnDown);
   bind_shortcut(fig, 'ctrl+q', @onDone);
   bind_shortcut(fig, 'f1', {@show_help, 'clustering', 'press ctrl+q when done'});
   
   % set(fig, 'closeRequestFcn', @onDone);
   
   % indicate preset channels
   for ch=make_vector(presetChannels, 'row')
      drawSelectionIndicator(hax, ch);
   end

   setappdata(fig, 'channelsToSubcluster', presetChannels);
   
   set(fig, 'UserData', 0);
   waitfor(fig,'UserData');
   channelsToSubcluster=getappdata(fig, 'channelsToSubcluster');
   
function onDone(sender, event) 
   fig=sender;
   set(fig, 'UserData', 1);   
   
function onAxesBtnDown(sender,eventdata)
   hax=sender;
   fig=get_parent_fig(hax);
   pos=get(hax, 'currentpoint');
   ch=round(pos(1));

   channelsToSubcluster=getappdata(fig, 'channelsToSubcluster'); 
   if(any(channelsToSubcluster==ch))
      return;
   end % if, already there
   
   drawSelectionIndicator(hax, ch);
   
   % NOTE: here the "channelsToSubcluster" are channel indices, not channel
   % numbers
   channelsToSubcluster(end+1)=ch;
   setappdata(fig, 'channelsToSubcluster', channelsToSubcluster);
   
   
   
function drawSelectionIndicator(hax, ch)
   x=ch;
   hLine=line([x x], get(hax, 'ylim'));
   set(hLine, 'LineStyle', ':', 'color', 'black');
   set(hLine, 'ButtonDownFcn', {@onLineBtnDown, struct('channel', ch)});

   
function onLineBtnDown(sender, eventdata, args)
   hLine=sender;
   ch=args.channel;
   fig=get_parent_fig(hLine);
   % remove the channel from the "result"
   channelsToSubcluster=getappdata(fig, 'channelsToSubcluster'); 
   channelsToSubcluster=setdiff(channelsToSubcluster, ch);
   setappdata(fig, 'channelsToSubcluster', channelsToSubcluster);
   
   % delete the gui obj
   delete(hLine);
   
   
function options=fillOptions(options)
   default_options('whatToPca', 'aboveThresh');
   options = default(options,'window_suppress_on',0);
   options = default(options,'win_locations',NaN);
   options = default(options,'abs_pathname_correct_on',0);
   options = default(options,'useMedianToClean',true);
   

   