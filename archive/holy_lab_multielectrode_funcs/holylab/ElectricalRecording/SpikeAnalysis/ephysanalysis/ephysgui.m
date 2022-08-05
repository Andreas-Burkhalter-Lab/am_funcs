function ephysgui(varargin)
% EPHYSGUI: Provide a context menu for ephys plots
% Calling syntaxes:
%   ephysgui(ephysin,plotparams)
%   ephysgui(axhandle,plotparams)
%   ephysgui(fighandle,plotparams)
%   ephysgui(fieldtoplot)
%
% Only the first of these is likely to be used by a user; the rest are
% for internal use.  Note that ephysgui differs from ephysplot in that
% any needed data will be fetched from disk for you.
%
% The input parameters:
%   ephys contains the input ephys data;
%   plotparams is a structure which specifies the initial plot type (see
%     EPHYSPLOTPARAMS);
%
% See also: EPHYS, EPHYSPLOT, EPHYSPLOTPARAMS, EPHYSFETCH.

% Copyright 2001 Timothy E. Holy, Jason Guo
% Changelog:
%   2004-05-03 (TEH): Fix channelnumber/cellnumber issue. I see Jason
%     tackled this already at the level of ephysgui, but this fix makes
%     use of the changes I'd made to ephysplot. I employed his idea of using
%     a dialog to set cellnumber and channelnumber.
%   2004-07-21 (TEH): For cell numbers, replace 'allonchan' with real
%     cell numbers early on. Also, the selection dialogs for channels and
%     cells should default to the current value.
%   2006-03-26 (RCH): added option that prevents ephysgui from
%     trying to load in additional information if a field is empty 
%     (plotparams.no_load = 1; default =0)

plottingcmds = {'stimulus','sniptimes','snippets','celltimes','envelope'};

% Create new plot
if (nargin > 1 & isstruct(varargin{1}))
  hax = gca;
  ephys = varargin{1};
  plotparams = varargin{2};
  if ~isfield(plotparams,'no_load')
	  plotparams.no_load = 0;
  end
  % Check to see if we need to load any data to get started
  if isfield(ephys,'sortfile')
    if ~isfield(ephys,'cellnums') || ~isnumeric(ephys(1).cellnums)
      ephys = ephysfetch(ephys,'cellnums');  % see TEH 2004-07-21
    end
  end
  loadindex = ephystofetch(ephys,plotparams.fieldtoplot);
  if (~isempty(loadindex) & (plotparams.no_load~=1))
    if isfield(ephys,plotparams.fieldtoplot)
      ephys(loadindex) = ephysfetch(ephys(loadindex),plotparams.fieldtoplot);
    else
      if (length(loadindex) ~= length(ephys))
        error('We have a weird problem...');
      end
      ephys = ephysfetch(ephys(loadindex),plotparams.fieldtoplot);
    end
  end
  % Store the data in the figure
  if isempty(get(gcf,'UserData'))
    set(gcf,'UserData',ephys);
  else
    warning('UserData wasn''t empty')
  end
  set(hax,'UserData',plotparams,'NextPlot','replacechildren');
  % Define the context menu
  cmenu = uicontextmenu;
  % Associate axes with the context menu
  set(hax,'UIContextMenu', cmenu);
  % Define callbacks for context menu items, depending
  % on what data is available or can be fetched
  datafiles = {'stimulusfile','snipfile','snipfile','sortfile','envelopefile'};
  submenu = {{'overlay','first','repeats','bar'},{'raster','PSTH','PSTH w/ sem'},{'overlay','reconstruct'},{'raster','PSTH','PSTH w/ sem'},{}};
  for i = 1:length(plottingcmds)
    if (isfield(ephys,plottingcmds{i}) | isfield(ephys,datafiles{i}))
      if isempty(submenu{i})
        cb = ['ephysgui(''',plottingcmds{i},''')'];
        item = uimenu(cmenu,'Label',plottingcmds{i},'Callback',cb);
      else
        item = uimenu(cmenu,'Label',plottingcmds{i});
        for j = 1:length(submenu{i})
          cb = ['ephysgui({''',plottingcmds{i},''',''',submenu{i}{j},'''})'];
          uimenu(item,'Label',submenu{i}{j},'Callback',cb);
        end
      end
    end
  end
  % Determine whether any "channel" data types are present
  channelfields = {'snipfile','envelopefile','wavefile','sniptimes','envelope','wave'};
  channelpresent = 0;
  for i = 1:length(channelfields)
    if isfield(ephys,channelfields{i})
      channelpresent = 1;
    end
  end
  if ~isfield(ephys,'channels')
      channelpresent = 0; % override in this case
  end
  % Determine whether any "cell" data types are present
  cellfields = {'sortfile','celltimes'};
  cellpresent = 0;
  for i = 1:length(cellfields)
    if isfield(ephys,cellfields{i})
      cellpresent = 1;
    end
  end
  % Set up the callback to call CASS
  if ((isfield(ephys,'sortfile') && isfield(ephys,'sort_cass') && ephys(1).sort_cass) ||...
      (isfield(ephys,'sort_template') && ephys(1).sort_template))
    uimenu(cmenu,'Label','View sorting',...
      'Callback',{@epgui_callcass,hax});
  end
  % Set up the remaining menu items
  defcmds = {'edit tags'};
  if supportlegacyppnumber
    defcmds{end+1} = 'number';
  else
    if channelpresent && length(unique([ephys.channels])) > 1
      defcmds{end+1} = 'channelnumber';
    end
    if cellpresent && length(unique([ephys.cellnums])) > 1
      defcmds{end+1} = 'cellnumber';
    end
  end
  defcmds = {defcmds{:},'slider window','manual'};
  for i = 1:length(defcmds)
    cb = ['ephysgui(''',defcmds{i},''')'];
    item = uimenu(cmenu,'Label',defcmds{i},'Callback',cb);
  end
  ephysplot(ephys,plotparams);
  return
end


% Refresh with a new plotparams structure
% Calling syntax: ephysgui(fighandle,plotparams)
%             or: ephysgui(axhandle,plotparams)
if (nargin > 1 & ishandle(varargin{1}))
  h = varargin{1};
  htype = get(h,'Type');
  if strcmp(htype,'figure')
    hax = findobj(varargin{1},'Tag','ephysplotax');
  elseif strcmp(htype,'axes')
    hax = varargin{1};
  end
  nax = length(hax);
  if (nax ~= length(varargin{2}))
    error('Number of axes is not equal to number of plotparams');
  end
  for i = 1:nax
    plotparams = varargin{2}(i);
    set(hax(i),'UserData',plotparams);
    axes(hax(i));
    ephysgui(plotparams.fieldtoplot);
  end
  return;
end

% Context menu callback or plotting command: must check to see if
% all the data we need are loaded, then can plot
if (nargin == 1 & (ischar(varargin{1}) | iscell(varargin{1})))
  hax = gca;
  plotparams = get(hax,'UserData');
  hfig = get(hax,'Parent');
  ephys = get(hfig,'UserData');
  isplottingcmd = 1;
  if ischar(varargin{1})
    if any(strmatch(varargin{1},plottingcmds,'exact'))
      plotparams.fieldtoplot = varargin{1};
    else
      isplottingcmd = 0;
    end
  else
    plotparams.fieldtoplot = varargin{1}{1};
    plotparams.type = varargin{1}{2};
  end
  if isplottingcmd
    % Check to see that the data has been loaded,
    % and load it if not
    % Find the ones to be plotted, and just load that data
    if ischar(plotparams.tags)
      tagindex = {strmatch(plotparams.tags,{ephys.tag},'exact')'};
      ntags = 1;
    else
      ntags = length(plotparams.tags);
      for i = 1:ntags
        tagindex{i} = strmatch(plotparams.tags{i},{ephys.tag},'exact')';
      end
    end
    alltags = [tagindex{:}];
    loadindex = alltags(ephystofetch(ephys(alltags),plotparams.fieldtoplot));
    if ~isempty(loadindex) & plotparams.no_load ~= 1
      ephystemp = ephysfetch(ephys(loadindex),plotparams.fieldtoplot);
      % Figure out which new fields have been created, and create space
      % for them in the proper order, so the subscripting command will work
      % (Note: some fetch routines save more data than just the new field name)
      [newfield,nford] = setdiff(fieldnames(ephystemp),fieldnames(ephys));
      [nfords,nfordi] = sort(nford);
      newfield = newfield(nfordi);
      for i = 1:length(newfield)
        eval(['[ephys.',newfield{i},'] = deal([]);']);   % Create the field
      end
      ephys(loadindex) = ephystemp;
      set(hfig,'UserData',ephys);  % Save new data in parent figure
    end
    % We may have to fill in channel plotting parameters
    pptmp = plotparams; % don't save channelnumber permanently
    if ~isempty(strmatch(plotparams.fieldtoplot,{'sniptimes','envelope','wave'},'exact')) && ...
        ~isfield(plotparams,'channelnumber')
      cellchandef = ephys(1).cellchandef{find(plotparams.cellnumber == ephys(1).cellnums)};
      pptmp.channelnumber = cellchandef(1);
    end
    ephysplot(ephys,pptmp);
  else
    % Not a plotting command, what is it?
    command = varargin{1};
    switch command
      case 'edit tags',
        if isfield(plotparams,'alltags')
          utags = plotparams.alltags;
        else
          utags = unique({ephys.tag});
          plotparams.alltags = utags;
        end
        [ctags,tagindx] = intersect(utags,plotparams.tags);
        showtags = utags;
        if isfield(plotparams,'tagstext')
          showtags = plotparams.tagstext;
        end
        [newtagindx,ok] = ...
          listdlg('ListString',showtags,'SelectionMode','multiple','InitialValue',tagindx);
        if ok
          plotparams.tags = utags(newtagindx);
          set(gca,'UserData',plotparams);
          ephysgui(plotparams.fieldtoplot);
          return;
        end
      case 'number',
        answer = inputdlg('Enter new cell number:');
        plotparams.number = str2num(answer{1});
        set(gca,'UserData',plotparams);
        ephysgui(plotparams.fieldtoplot);
        return;
      case 'channelnumber',
        ttChannels=ephys(1).channels;
        ttStrs=num2str(sort(ttChannels'));
        [sel,ok]=listdlg('promptstring', 'please select a channel', ...
                         'selectionmode', 'single', ...
                         'name', 'ephysgui', ...
                         'liststring', ttStrs, ...
                         'InitialValue',find(sort(ephys(1).channels) == plotparams.channelnumber) ...
                         );
        if(ok==0) return; end
        ttCh=str2num(ttStrs(sel,:));
        plotparams.channelnumber = ttCh;
        set(gca,'UserData',plotparams);
        ephysgui(plotparams.fieldtoplot);
        return;
      case 'cellnumber',
        ttCells=ephys(1).cellnums;
        ttStrs=num2str(sort(ttCells'));
        [sel,ok]=listdlg('promptstring', 'please select a cell', ...
                         'selectionmode', 'single', ...
                         'name', 'ephysgui', ...
                         'liststring', ttStrs, ...
                         'InitialValue',find(sort(ephys(1).cellnums) == plotparams.cellnumber) ...
                         );
        if(ok==0) return; end
        ttCl=str2num(ttStrs(sel,:));
        plotparams.cellnumber = ttCl;
        set(gca,'UserData',plotparams);
        ephysgui(plotparams.fieldtoplot);
        return;
      case 'slider window',
        swoptions.restoreonclose = 1;
        pos = get(gcf,'Position');
        sliderwindow(gca,swoptions); % Old position: [16 50 988 pos(4)]
      case 'manual',
        keyboard
      otherwise
        error('Not yet implemented');
    end
  end
else
  error('Unknown calling syntax');
end
set(gca,'UserData',plotparams); % Save the plotparams in the current axis

function epgui_callcass(hobject,eventdata,hax)
  hfig = get_parent_fig(hobject);
  ephys = get(hfig,'UserData');
  plotparams = get(hax,'UserData');
  if ~isfield(ephys,'sort_template')
    sortfile = ephys(1).sortfile;
    if (plotparams.cellnumber == round(plotparams.cellnumber))
      % The cass channel is the channel number
      cellchandef = ephys(1).cellchandef{find(plotparams.cellnumber == ephys(1).cellnums)};
      cass_chan = cellchandef(1);
    else
      % The cass channel is encoded in the cell number
      cass_chan = round(plotparams.cellnumber);
    end
  else
    % Special measures are required for template-based sorting, to handle
    % the different part directories and multi-channel cell definition
    cellIndex = find(ephys(1).cellnums == plotparams.cellnumber);
    celltag = ephys(1).celltags{cellIndex};
    [sortfile,R] = strtok(celltag,';');
    cass_chan = str2num(R(2:end))
    cass_chan = cass_chan(1);
  end
  cass('dirname',sortfile,'channel',cass_chan);
  
  