function varargout = clean_snipfiles(varargin)
% CLEAN_SNIPFILES: eliminate weird spikes from snipfiles
% Syntax:
%   clean_snipfiles(filenames)
% where
%   filenames is a cell array of .ssnp filenames.
%
% See also: AUTOSORT, CASS.

% Copyright 2005 by Timothy E. Holy

% CLEAN_SNIPFILES M-file for clean_snipfiles.fig
%      CLEAN_SNIPFILES, by itself, creates a new CLEAN_SNIPFILES or raises the existing
%      singleton*.
%
%      H = CLEAN_SNIPFILES returns the handle to a new CLEAN_SNIPFILES or the handle to
%      the existing singleton*.
%
%      CLEAN_SNIPFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLEAN_SNIPFILES.M with the given input arguments.
%
%      CLEAN_SNIPFILES('Property','Value',...) creates a new CLEAN_SNIPFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before clean_snipfiles_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to clean_snipfiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help clean_snipfiles

% Last Modified by GUIDE v2.5 12-Jul-2005 17:18:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clean_snipfiles_OpeningFcn, ...
                   'gui_OutputFcn',  @clean_snipfiles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before clean_snipfiles is made visible.
function clean_snipfiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to clean_snipfiles (see VARARGIN)

  % Choose default command line output for clean_snipfiles
  handles.output = hObject;

  if (nargin > 3 && length(varargin) > 0)
    filenames = varargin{1};
  else
    filenames = -1;
  end

  % See if a results file exists (this stores user choices)
  if (exist('clean_snip_temp_results') == 2)
    load clean_snip_temp_results
    if exist('landmarkKeep')
      % This makes sure that our file written to test permissions doesn't
      % end up fooling us
      clean_snipfiles_execute(sorthead,channels,landmarkKeep)
      clean_snipfiles_cleanup(hObject)
      return  % all done!
    end
  end

  % See if a temp sorting file exists
  use_temp_file = 0;
  if exist('clean_snip_temp','dir')
    % See if it contains the same files supplied on the command line
    if isequal(filenames,-1)
      use_temp_file = 1;
    else
      load('clean_snip_temp/overview')
      fh = [sorthead.fh];
      for i = 1:length(fh)
        temp_filenames{i} = fh(i).filename;
      end
      if isequal(filenames,temp_filenames)
        use_temp_file = 1;
      else
        warning(['Files in clean_snip_temp did not match supplied files, '...
                 'starting from scratch'])
      end
    end
  end

  % If no pre-existing temp file, do the calculations and save to temp file
  if ~use_temp_file
    sorthead = snipfile2sortheader(filenames);
    n_landmarks = 36;   % max 6-by-6 panel of axes
    options_autosort = struct('t2V',0,...
                              'cluster_func',@meanshift,...
                              'use_projection',0,...
                              'n_landmarks',n_landmarks,...
                              'max_snips_to_cluster',10000,...
                              'Rfactor',.2); % don't overcluster! (changed even lower - RCH) 
    autosort(sorthead,'clean_snip_temp',options_autosort);
  end
  
  % Load data from temp file
  load('clean_snip_temp/overview')
  channels = sorted_get_channels('clean_snip_temp');
  
  % Draw grid of axes
  subplotdims = ceil(sqrt(options_autosort.n_landmarks));
  subplotdims = subplotdims([1 1]);
  base_position = get(handles.axesShowWaveforms,'Position');
  widthTotal = 0.9/subplotdims(2);
  heightTotal = 0.9/subplotdims(1);
  leftmargin = 0.05;
  bottommargin = 0.05;
  widthfrac = 0.95;  % affects the gap size between panels
  heightfrac = 0.95;
  for i = 1:options_autosort.n_landmarks
    [row,col] = ind2sub(subplotdims,i);
    row = subplotdims(1)-row+1;
    norm_position = [leftmargin+(col-1)*widthTotal ...
                     bottommargin+(row-1)*heightTotal ...
                     widthfrac*widthTotal ...
                     heightfrac*heightTotal];
    position = SetAxPosNorm(base_position,norm_position);
    handles.panel(i) = axes('Parent',handles.figCleanSnipFiles,...
                            'Position',position,...
                            'XTick',[],...
                            'YTick',[],...
                            'SelectionHighlight','on');
    install_mouse_event_handler(handles.panel(i),'up',@axes_mouse_up);
  end
  install_mouse_event_handler(handles.panel(i),'up',@figure_mouse_up);

  % Store data in figure
  setappdata(handles.figCleanSnipFiles,'sorthead',sorthead);
  setappdata(handles.figCleanSnipFiles,'channels',channels);
  setappdata(handles.figCleanSnipFiles,'channelIndex',1);
  setappdata(handles.figCleanSnipFiles,'options_autosort',options_autosort);

  % Update handles structure
  guidata(hObject, handles);

  clean_snipfiles_load_channel(handles)
  clean_snipfiles_display_channel(handles)
    

% UIWAIT makes clean_snipfiles wait for user response (see UIRESUME)
% uiwait(handles.figCleanSnipFiles);


% --- Outputs from this function are returned to the command line.
function varargout = clean_snipfiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if (isfield(handles,'output') && ishandle(handles.output))
  varargout{1} = handles.output;
end


function result = axes_mouse_up(sender, event_args)
  is_shift_down=strcmp(get(get_parent_fig(sender), 'SelectionType'), ...
                       'extend');
  if(is_shift_down) % toggle current axes only
    if(strcmp(get(sender,'selected'), 'on'))
      set(sender, 'selected', 'off');
    else
      set(sender, 'selected', 'on');
    end
  else
    handles=guidata(sender);
    set(handles.panel, 'selected', 'off');
    set(sender, 'selected', 'on');
  end
  result = 1;
 

function result = figure_mouse_up(sender, event_args)
  handles = guidata(sender);
  set(handles.panel, 'selected', 'off');
  result = 1;
  

function clean_snipfiles_load_channel(handles)
  sorthead = getappdata(handles.figCleanSnipFiles,'sorthead');
  channels = getappdata(handles.figCleanSnipFiles,'channels');
  channelIndex = getappdata(handles.figCleanSnipFiles,'channelIndex');
  options_autosort = getappdata(handles.figCleanSnipFiles, ...
                                'options_autosort');
  dirChan = ['clean_snip_temp' filesep 'chan' ... 
             num2str(channels(channelIndex)) filesep];
  sortchan = load_sort_info([dirChan 'autosort_info']);
  % Store the results
  setappdata(handles.figCleanSnipFiles,'landmarkPosition', ...
                           sortchan.landmarkWaveform);
  setappdata(handles.figCleanSnipFiles,'landmarkClust', ...
                           sortchan.landmarkClust);
  % Update display counter
  set(handles.textChannel,'String',...
    sprintf('Channel %d (%d/%d)',channels(channelIndex),channelIndex,...
    length(channels)));


function clean_snipfiles_display_channel(handles)
  landmarkPosition = getappdata(handles.figCleanSnipFiles,'landmarkPosition');
  landmarkClust = getappdata(handles.figCleanSnipFiles,'landmarkClust');
  sorthead = getappdata(handles.figCleanSnipFiles,'sorthead');
  % Display landmark waveforms in panels.
  uclust = setdiff(unique(landmarkClust),0);  % Skip the trash
  nclust = length(uclust);
  clabel = cell(1,nclust);
  for i = 1:nclust
    clabel{i} = find(landmarkClust == uclust(i));
  end
  sniprange = unique([sorthead.sniprange]);
  x = sniprange(1):sniprange(2);
  allShownLandmarksIndex = cat(2,clabel{:});
  allShownLandmarks = landmarkPosition(:,allShownLandmarksIndex);
  ylim = [min(min(allShownLandmarks)) max(max(allShownLandmarks(:)))];
  for i = 1:nclust
    plot(handles.panel(i),x,landmarkPosition(:,clabel{i}));
    setappdata(handles.panel(i),'landmarkIndex',clabel{i});
    % Plot line at peak
    line([0 0],ylim,'Parent',handles.panel(i),...
      'Color','k','LineStyle','--');
  end
  % Clear the remaining axes
  for i = nclust+1:length(handles.panel)
    cla(handles.panel(i))
    if isappdata(handles.panel(i),'landmarkIndex');
      rmappdata(handles.panel(i),'landmarkIndex');
    end
  end
  set(handles.panel,'XTick',[],'YTick',[],...
    'XLim',x([1 end]),'YLim',ylim)
  

function clean_snipfiles_execute(sorthead,channels,landmarkKeep)
  % We're not assuming the GUI exists anymore, because we might be
  % parsing a temp file
  nchannels = length(channels);
  for i = 1:length(sorthead)  % Do one file at a time
    filein = sorthead(i).fh.filename;
    [pathstr,basename,ext] = fileparts(filein);
    fileout = [basename '_clean' ext];
    % Test for write permission errors now before we have invested a lot
    % of time; at the same time, take the opportunity to write the header
    % to the new snippet file
    strHeader = read_ascii_header(filein);  % read old header
    strHeader = update_value(strHeader,'channel list',num2str(channels));
    dummystr = num2str(zeros(size(channels)));
    strHeader = update_value(strHeader,'timesfpos',dummystr);
    strHeader = update_value(strHeader,'snipsfpos',dummystr);
    strHeader = update_value(strHeader,'finetimesfpos',dummystr);
    strHeader = update_value(strHeader,'detpeaksfpos',dummystr);
    [fidout,msg] = fopen(fileout,'w');
    if (fidout < 0)
      error(['Can''t open file ' fileout ' for writing; check ' ...
                          'permissions?']);
    end
    update_header(fidout,strHeader);
    for j = 1:nchannels
      progstr = sprintf('Saving results on file %s (%d/%d), channel %d (%d/%d)',...
                        filein,i,length(sorthead),...
                        channels(j),j,nchannels);
      progress_bar(struct('progress',(i-1)*nchannels+j,...
                          'max',length(sorthead)*nchannels,...
                          'what',progstr));
      shc = sortheader_importchan(sorthead(i),channels(j));
      nsnips = shc.numofsnips;
      % Process snippets in blocks to avoid overfilling memory
      blocksize = 5000;
      nblocks = ceil(nsnips/blocksize);
      tKeepflag = cell(1,nblocks);
      for k = 1:nblocks
        % Read in a block of snippets
        startIndex = (k-1)*blocksize+1;
        endIndex = min(nsnips,startIndex+blocksize-1);
        snip = sortheader_readsnips(shc,startIndex:endIndex);
        % Determine the closest landmark
        [d,landmarkIndex] = mindist(snip,landmarkKeep(j).position);
        % Keep the snippets that are closest to the valid landmarks
        tKeepflag{k} = landmarkKeep(j).keepflag(landmarkIndex);
      end
      % Concatenate blocks together
      keepflag = cat(2,tKeepflag{:});
      % Set up the mmap structure for convenient non-overflowing output
      snipmm = snipfile_mmap(filein,channels(j));
      % Write the snippets to a file
      strHeader = snipfile_append_channel(fidout,strHeader,channels(j),snipmm,keepflag);
    end
  end

function clean_snipfiles_cleanup(fighandle)
  delete('clean_snip_temp_results.mat')
  rmdir('clean_snip_temp','s')
  if ishandle(fighandle)
    close(fighandle)
  end
  

% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% hObject    handle to btnNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Store the results
  if isappdata(handles.figCleanSnipFiles,'landmarkKeep')
     landmarkKeep = getappdata(handles.figCleanSnipFiles,'landmarkKeep');
   else
     landmarkKeep = [];
   end
  landmarkPosition = getappdata(handles.figCleanSnipFiles,'landmarkPosition');
  landmarkClust = getappdata(handles.figCleanSnipFiles,'landmarkClust');
  landmarkKeep(end+1).position = landmarkPosition;
  landmarkKeep(end).keepflag = (landmarkClust > 0);
  setappdata(handles.figCleanSnipFiles,'landmarkKeep',landmarkKeep);
  % Advance the channel index
  channels = getappdata(handles.figCleanSnipFiles,'channels');
  channelIndex = getappdata(handles.figCleanSnipFiles,'channelIndex');
  channelIndex = channelIndex+1;
  setappdata(handles.figCleanSnipFiles,'channelIndex',channelIndex);
  if (channelIndex > length(channels))
    % We're done! Save the results for safety
    sorthead = getappdata(handles.figCleanSnipFiles,'sorthead');
    try
      save clean_snip_temp_results sorthead channels landmarkKeep
    catch
      warning('Unable to save results file')
    end
    % Now implement the decisions
    clean_snipfiles_execute(sorthead,channels,landmarkKeep)
    % Delete the temp files
    clean_snipfiles_cleanup(handles.figCleanSnipFiles)
  elseif (channelIndex == length(channels))
    set(handles.btnNext,'String','Done')
  end

  
% --- Executes on key press over figCleanSnipFiles with no controls selected.
function figCleanSnipFiles_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figCleanSnipFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  currentChar = get(hObject,'CurrentCharacter');
  if strcmp(currentChar,char(8))
    selected = get(handles.panel,'Selected');
    selectedIndex = strmatch('on',selected,'exact');
    if ~isempty(selectedIndex)
      landmarkClust = getappdata(handles.figCleanSnipFiles,'landmarkClust');
      for i = 1:length(selectedIndex)
        landmarkClust(getappdata(handles.panel(selectedIndex(i)),'landmarkIndex')) = 0;
      end
      setappdata(handles.figCleanSnipFiles,'landmarkClust', ...
                               landmarkClust);
      set(handles.panel, 'selected', 'off');
      clean_snipfiles_display_channel(handles)
    end
  end

