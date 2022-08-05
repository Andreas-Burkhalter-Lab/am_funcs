function deltarategui(varargin)
   if nargin > 1
      % Initial setup
      ephys = varargin{1};
      if isfield(ephys,'celltimes')
        fieldtoplot = 'celltimes';
        ppidname = 'cellnumber';
      else
        fieldtoplot = 'sniptimes';
        ppidname = 'channelnumber';
      end
      dr = varargin{2};
      drerr = varargin{3};
      tags = varargin{4};
      settingcolors = 0;
      usingtagstext = 0;
      itemsmarked = 0;
      if (nargin > 4)
         nextarg = 5;
         if isnumeric(varargin{5})
            col = varargin{5};
            settingcolors = 1;
            nextarg = 6;
         end
         if (nargin >= nextarg & iscell(varargin{nextarg}))
            tagstext = varargin{nextarg};
            usingtagstext = 1;
            nextarg = nextarg+1;
         end
         if (nargin >= nextarg)
            markers = varargin{nextarg};
            itemsmarked = 1;
         end
      end
      % Check special syntax for single files
      if ~iscell(dr)
         dr = {dr};
      end
      if ~iscell(ephys)
         ephys = {ephys};
      end
      nfiles = length(dr);
      ntags = length(tags);
      if (length(ephys) ~= nfiles)
         error('File number mismatch');
      end
      % Create zoom figures for each file in ephys, and set the UserData
      hfigs = zeros(1,nfiles);
      stimparams = struct('fieldtoplot','stimulus','type','overlay', ...
                          'channelnumber',1,'showtags',1,'axtype','stimulus','alltags',{tags},...
                          'axisproperties',{{'XTick',[]}},'titlenumber',1);
      dataparams = struct('fieldtoplot',fieldtoplot,'type','PSTH w/ sem', ...
                          'tomicrovolts',68,'showtags',1,'axtype','data',ppidname,1,...
                          'alltags',{tags});
      if (settingcolors)
         stimparams.tagcolors = col;
         dataparams.tagcolors = col;
      end
      if (usingtagstext)
         stimparams.tagstext = tagstext;
         dataparams.tagstext = tagstext;
      end
      for i = 1:nfiles
         hfigs(i) = figure('Name',ephys{i}(1).basefilename, ...
                           'NumberTitle','off','Visible','off',...
                           'CloseRequestFcn','deltarategui(''closezoomfig'')');
         % Make a stimulus axis and a data axis
         hax = SplitVert([0.85 0.88]); delete(hax(2)); hax = hax([1 3]);
         %set(hax(1),'UserData',stimparams);
         %set(hax(2),'UserData',dataparams);
         axes(hax(1));
         ephysgui(ephys{i},stimparams);
         axes(hax(2));
         ephysgui(ephys{i},dataparams);
         set(hfigs(i),'Visible','off')
      end
      
      % Create the user data cell arrays for marked items
      if itemsmarked
         markedfields = fieldnames(markers);
         ud = cell(1,nfiles);
         for i = 1:nfiles
            nchans = size(dr{i},2);
            ud{i} = cell(1,nchans);
            for j = 1:nchans
               ud{i}{j} = {};
            end
            for j = 1:length(markedfields)
               curmarked = markers(i).(markedfields{j});
               for k = 1:length(curmarked)
                  [com,indx] = intersect(ephys{i}(1).channels,curmarked(k));
                  if isempty(com)
                     error('Marked channel not recorded');
                  end
                  ud{i}{indx}{end+1} = markedfields{j};
               end
            end
         end
      end
      
      % Create the deltarate figure
      figure('CloseRequestFcn','deltarategui(''cleanup'')');
      set(gca,'Position',[0.21 .11 .775 .815]);
      breaks = (1:ntags+itemsmarked-1)/(ntags+itemsmarked);
      breaks = sort([breaks breaks+0.1*breaks(1)]);
      %breaks = breaks(1:end-1);
      haxdrf = SplitVert(breaks);
      delete(haxdrf(2:2:end));
      haxdrf = haxdrf(1:2:end);
      for i = 1:ntags
         %subplot(ntags+itemsmarked,1,i+itemsmarked);
         axes(haxdrf(i));
         hold on
         curx = 1;
         dividers = zeros(1,nfiles);
         for j = 1:nfiles
            [ntgs,nchan] = size(dr{j});
            for k = 1:nchan
               if(nfiles==1) % for Francesco's dodeltar.m, in that case nfiles is always 1.
                  % @note: here comes the limitation:
                  %    we assume all files specified in ephys have same channels field.
                  %    So ephys{1}(1) is used to convert channel idx to channel.
                  if strcmp(ppidname,'cellnumber')
                    ttCh = ephys{1}(1).cellnums(k);
                  else
                    ttCh=ephys{1}(1).channels(k);
                  end
                  if(str2num(version('-release'))>=14)
                    hline = errorbar('v6', ttCh,dr{j}(i,k),drerr{j}(i,k));  
                  else
                    hline = errorbar(ttCh,dr{j}(i,k),drerr{j}(i,k));
                  end
               else
                  if(str2num(version('-release'))>=14)
                    hline = errorbar('v6', curx,dr{j}(i,k),drerr{j}(i,k));  
                  else
                    hline = errorbar(curx,dr{j}(i,k),drerr{j}(i,k));
                  end
               end
               %set(hline(2),'Marker','o','MarkerSize',9);
               if settingcolors
                  set(hline,'Color',col(i,:));
               end
               cbstruct = struct('tag',tags{i},'fileindx',j, ...
                                 ppidname,ttCh);  % TEH 2004-05-03
               set(hline,'ButtonDownFcn','deltarategui(''focus'')',...
                         'UserData',cbstruct);
               if itemsmarked
                  tmp = ud{j}{k};
                  if ~isempty(tmp)
                     [tmp{2,:}] = deal(':');
                     set(hline,'Tag',[':',[tmp{:}]]);
                  else
                     set(hline,'Tag','');
                  end
               end % if,
               curx = curx+1;
            end % for, k = 1:nchan 
            dividers(j) = curx;
            curx = curx+1;
         end % for, j = 1:nfiles
         hold off
         axis tight
         set(gca,'TickDir','out','XLim',[0 curx]);
         ylim = get(gca,'YLim');
         hdline = line(repmat(dividers(1:end-1),2,1),repmat(ylim',1,nfiles-1));
         set(hdline,'Color',[0 0 0],'LineStyle',':');
         line([0 curx],[0 0],'Color',[0 0 0],'LineStyle','--');
         if (i < ntags)
            set(gca,'XTick',[]);
         end % if,
         %title(tags{i});
         hyl = ylabel(tags{i},'Rotation',0,'HorizontalAlignment','Right',...
           'VerticalAlignment','Middle');
      end % for i = 1:ntags
      figdata = struct('hfigs',hfigs);
      figdata.ppidname = ppidname;
      if itemsmarked
         % Save the name of the markers
         figdata.markernames = markedfields;
         % Draw button to select marked items
         subplot(ntags,1,1);
         pos = get(gca,'position');
         pos(2) = pos(2) + pos(4)/2;
         pos(4) = pos(4)/2;
         delete(gca);
         hb1 = uicontrol(gcf,'Style','PushButton','String','Select Markers',...
                         'units','normalized','Position',pos,...
                         'Callback','deltarategui markers');    
      end
      set(gcf,'UserData',figdata);  % Store the figure handles for later use
                                    % Set up the button for 
   else
      switch varargin{1},
         case 'focus',
            % Determine which object was clicked, and which figure to open
            cbstruct = get(gcbo,'UserData');
            figdata = get(gcbf,'UserData');
            hfigs = figdata.hfigs;
            hfig = hfigs(cbstruct.fileindx);
            % Get the axes to plot into
            hax = findobj(hfig,'Type','axes');
            ud = get(hax,'UserData');
            for i = 1:length(ud)
               if (isstruct(ud{i}) & isfield(ud{i},'axtype'))
                  ud{i}.tags = cbstruct.tag;
                  ud{i}.(figdata.ppidname) = cbstruct.(figdata.ppidname);  % TEH 2004-05-03, 2007-03-17
                  ephysgui(hax(i),ud{i});
               end
            end
            set(hfig,'Visible','on');
         case 'closezoomfig',
            selection = questdlg('Really close figure?',...
                                 'Close Request Function',...
                                 'Yes','No','Yes');
            switch selection,
               case 'Yes',
                  delete(get(0,'CurrentFigure'))
               case 'No'
                  return
            end
         case 'markers'
            figdata = get(gcbf,'UserData');
            [newmarkers,ok] = listdlg('ListString',figdata.markernames,...
                                      'SelectionMode','multiple','InitialValue',[]);
            if ok
               alllines = findobj(gcbf,'Type','line');
               set(alllines,'LineWidth',1);
               linetags = get(alllines,'Tag');
               for i = 1:length(newmarkers)
                  indx = [];
                  for j = 1:length(alllines)
                     if findstr([':',figdata.markernames{newmarkers(i)},':'],linetags{j})
                        indx(end+1) = j;
                     end
                  end
                  set(alllines(indx),'LineWidth',4);
               end
               %         hline = findobj(gcbf,'Tag','deltarate');
               %         for i = 1:length(hline)
               %           clinetag = get(hline(i),'UserData');
               %           clinetag = clinetag.marks;
               %           if ~isempty(intersect(clinetag,figdata.markernames(newmarkers)))
               %             set(hline,'LineWidth',2);
               %           else
               %             set(hline,'LineWidth',1);
               %           end
               %         end
            end
         case 'cleanup'
            thisfig = get(0,'CurrentFigure');
            ud = get(thisfig,'UserData');
            if isstruct(ud)
              hfigs = ud.hfigs;
              hfigs = hfigs(find(ishandle(hfigs)));
              set(hfigs,'CloseRequestFcn','delete(get(0,''CurrentFigure''))');
            end
            delete(thisfig);
         otherwise
            error('Unrecognized calling syntax');
      end
   end
   
