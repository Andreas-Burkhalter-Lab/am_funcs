function plot_valve_transition(roiIntensitiesAxes)
   % get valve interested
   mouseDownPos=get(roiIntensitiesAxes, 'currentpoint');
   mouseDownVlvValue=mouseDownPos(1);
   valveToPlotTransition=round(mouseDownVlvValue);

   % appdata needed
   roi_label=getappdata(roiIntensitiesAxes, 'roi_label');
   figSummary=getappdata(roiIntensitiesAxes, 'figSummary');
   figureRois=getappdata(roiIntensitiesAxes, 'figureRois');

   ip=getappdata(figSummary, 'ip');

   if(~isfield(ip, 'raw'))
      uiwait(msgbox('You must provide argument raw in imstimviewdf() to make ploting valve transition possible.','Manage ROIs','modal'));
      return;
   end

   % find out the related valve figure:
   valveFigures=getappdata(figSummary, 'valveFigures');
   validVlvFig=[]; validValves=[];
   for idxValveFig=1:length(valveFigures)
      ttFigHandle=valveFigures(idxValveFig);
      if(ishandle(ttFigHandle) && ~isempty(getappdata(ttFigHandle, 'valve')))
         validVlvFig(end+1)=ttFigHandle;
         validValves(end+1)=getappdata(validVlvFig(end), 'valve');
      end
   end
   tIdx=find(validValves==valveToPlotTransition);
   if(isempty(tIdx)) return; end
   relatedValveFig=validVlvFig(tIdx);

   allAxesInVlvFig=getappdata(relatedValveFig, 'allAxes');
   
   % Loop over trials
   for idxTrial=1:length(allAxesInVlvFig)
     progress_bar(struct('progress',idxTrial,'max',length(allAxesInVlvFig),'what', ['Trial #' num2str(idxTrial)]));
     % Determine which frames should be analyzed
     axesEachTrial=allAxesInVlvFig(idxTrial);
     stacknum = getappdata(axesEachTrial, 'movie_stacknum');
     movie_filename = getappdata(axesEachTrial, 'movie_filename');
     fileIndex = strmatch(movie_filename,{ip.raw.filename},'exact');
     movie_indices{idxTrial}=fileIndex(find([ip.raw(fileIndex).stacknum] >= stacknum(1) & ...
                                  [ip.raw(fileIndex).stacknum] <= stacknum(2)));
     trialNums{idxTrial}=getappdata(axesEachTrial, 'trial');
     movie_indices_len(idxTrial)=length(movie_indices{idxTrial});
     tIp = ip.raw(movie_indices{idxTrial});
     % Get ROI definitions (this allows for drift, since ROI is
     % obtained at time just before valve transition), and pick the one we
     % want
     rois = get_roi(figSummary,stacknum(1));
     cellArrayRoi = split_roi_defs(rois);
     indxRoi = find(rois.label == roi_label);
     tRoi = cellArrayRoi{indxRoi};
     % Measure intensity
     intensity{idxTrial} = roimeasure({tIp.image},tRoi);
     stimulus{idxTrial} = [tIp.stimulus];
   end
   
   % Prepare to plot
   hFigTransition=figure('NumberTitle','off', ...
      'Name',['Valve ' num2str(valveToPlotTransition) ' transition details on roi ' num2str(roi_label)], ...
      'position', [100 100 500 700]);
%    [100 100 550 800]
%    axesRoiMeasurement=subplot(3,1,[1 2]);
%    axesStimuli=subplot(3,1,[3]);
%    set([axesRoiMeasurement axesStimuli], 'NextPlot', 'add');
      
   % intensities=reshape(intensities, length(movie_indices{1}), length(movie_indices)); % #cols=#trials
   % hLines=plot(intensites);
   % note: reshape() doesn't work b/c different trials may have different length
   tToAverage={};
   trialmap = trialcolor(length(movie_indices_len));
   hold on
   for idxTrial=1:length(movie_indices_len)
      tStim = stimulus{idxTrial};
      tTransitionIdx = [0 find(diff(tStim)) length(tStim)];
      bgIdx = 1:tTransitionIdx(2);
      tInten = intensity{idxTrial};
      bg = mean(tInten(bgIdx));
      tInten = (tInten-bg)/bg;  % deltaf/f
      %tColor=pick_next_color(axesRoiMeasurement);
      tColor = trialmap{idxTrial};
      %hLines(idxTrial)=plot(axesRoiMeasurement, tInten,'color',tColor,'linewidth',1.5);
      hLines(idxTrial)=plot(tInten,'color',tColor,'LineWidth',0.25);
      tToAverage{idxTrial}=tInten;
      %hStairs(idxTrial)=stairs(axesStimuli, tStim);
      % Plot stimulus as bar plot
%       for idxStim = 2:length(tTransitionIdx)
%         if tStim(tTransitionIdx(idxStim))            % if not flush
%           xr = [tTransitionIdx(idxStim-1)+1 tTransitionIdx(idxStim)];
%           yr = [-0.25 0.25]+idxTrial;
%           patch(xr([1 2 2 1]),yr([1 1 2 2]),tColor);
%         end
%       end
      %set([hLines(idxTrial) hStairs(idxTrial)], 'color', tColor);
%       Plot stimulus as bar plot
if idxTrial == 1
    for idxStim = 2:length(tTransitionIdx)
        if tStim(tTransitionIdx(idxStim))            % if not flush
            xr = [tTransitionIdx(idxStim-1)+1 tTransitionIdx(idxStim)];
        end
    end
end
      
   end
   
   hax = gca;
   ylim = get(hax, 'ylim');
   lowbound = ylim(1);
   yr = [lowbound lowbound*.95];
   patchcolor = [0 0 0];
   patch(xr([1 2 2 1]),yr([1 1 2 2]),patchcolor);

   %    trialNumsAsStr=cell(1, length(trialNums));
   %    [trialNumsAsStr{:}]=foreach_g(@num2str, trialNums{:});
%    legend(hLines, trialNumsAsStr, 'location', 'NorthWest');
   
   % plot the average:
   average=noneven_mean(tToAverage);
   % tColor=pick_next_color(axesRoiMeasurement);
   %hLineAverage=plot(axesRoiMeasurement, average, 'color', 'black');
   hLineAverage=plot(average, 'color', 'black');
   set(hLineAverage, 'linewidth', 3, 'linestyle', '-'); 
hold off

xlim = [1 max(movie_indices_len)];
set(gca,'xlim',xlim);

   % add context menu for re-calc mean
   cmenu = uicontextmenu;
   set(hax,'UIContextMenu',cmenu);
   uimenu(cmenu,'Label','Recalculate mean','Callback',{@fRecalcMean, struct('hLineAverage', hLineAverage)});
   uimenu(cmenu,'Label','Refilter ...','Callback',{@fRefilter, struct('hLineAverage', hLineAverage)});

   xlabel('Frames','FontSize',12,'FontWeight','bold');
   ylabel('DF/F (percent change)','FontSize',12,'FontWeight','bold');
%    xlabel(hFigTransition, 'Frames');
%    xlabel(axesStimuli, 'Frames');
   first_ipraw=ip.raw(movie_indices{1}(1) );
   if(isfield(first_ipraw, 'computation'))
%       ylabel(hFigTransition, first_ipraw.computation);
      ylabel(first_ipraw.computation);
   end
   %ylabel(axesStimuli, 'Stimuli');
   % ylabel(axesStimuli, 'Trial');
   xlim=get(hax, 'xlim');
   set(gca, 'xlim', [1 xlim(2)],'ylim',[ylim(1) ylim(2)],'TickDir','out','box','off')
%    set([axesRoiMeasurement axesStimuli], 'xlim', [1 xlim(2)]);
%    set(axesStimuli,'YDir','reverse','YTick',1:length(movie_indices_len))
   
function fRecalcMean(sender, event_data, args)
   % sender is the uimenu
   hLineAverage=args.hLineAverage;
   
   % find lines to average excluding hLineAverage and those who have * marker 
   tAxes=get(hLineAverage, 'parent');
   linesToAverage=findobj(tAxes, 'type', 'line');
   linesToExclude=findobj(linesToAverage, 'flat', 'Marker', '*');
   linesToAverage=setdiff(linesToAverage, [linesToExclude(:); hLineAverage]);
   
   % calc the average (the xdata calculation is a lazy way to get the
   % longest xdata. )
   valuesToAverage=get(linesToAverage, 'ydata');
   xdata_all=get(linesToAverage, 'xdata');
   if(~iscell(valuesToAverage))
      valuesToAverage={valuesToAverage};
      xdata_all={xdata_all};
   end
   average=noneven_mean(valuesToAverage);
   xdata=noneven_mean(xdata_all);
   
   % re-set the xdata and ydata of hLineAverage
   % xdata=get(hLineAverage, 'xdata');
   set(hLineAverage, 'xdata', xdata(1:length(average)), 'ydata', average);
   
function fRefilter(sender, event_data, args)
   % sender is the uimenu
   hLineAverage=args.hLineAverage;
   
   % find all lines excluding the average line
   tAxes=get(hLineAverage, 'parent');
   lines=findobj(tAxes, 'type', 'line');
   lines=setdiff(lines, hLineAverage);
   
   % get the N for medfilt1
   dlgTitle='Please input filter scope';
   prompt={'filter scope (the bigger, the smoother. Input 1 for no filtering.):', ...
           };
   defaultValues={num2str(1)};          
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,defaultValues);
   if(isempty(answer)) return ; end
   filterscope=str2num(answer{1});
   
   for idxLine=1:length(lines)
      curLine=lines(idxLine); 
      oldData=getappdata(curLine, 'oldData');
      if(isempty(oldData))
         oldData=get(curLine, 'ydata');
         setappdata(curLine, 'oldData', oldData);
      end
      set(curLine, 'ydata', medfilt1(oldData,filterscope));
   end
   
   % update the average also:
   fRecalcMean(sender, [], args);
