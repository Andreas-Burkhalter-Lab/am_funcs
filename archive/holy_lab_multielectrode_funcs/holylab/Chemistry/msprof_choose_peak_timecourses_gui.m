function varargout = msprof_choose_peak_timecourses_gui(varargin)
% MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI
% Syntax:
%    msprof_choose_peak_timecourses_gui
%    msprof_choose_peak_timecourses_gui(basename)
%  These forms start a new instance of the GUI.
%
%    msprof_choose_peak_timecourses_gui(mz)
%    msprof_choose_peak_timecourses_gui(mz,trange)
%  This syntax causes an existing instance of the GUI to "navigate" to a
%  given m/z (and if supplied, visually indicate a particular trange)
%
%    msprof_choose_peak_timecourses_gui(...,options)
%  With either of the above, you can add some options.
%
% In the "middle", between the temporal plot and the axis for the raw
% intensity matrix data, is a window that shows the timespans that identify
% individual compounds.  You can manually modify these timespans:
%  1. To add a new one, click-drag across the region in the temporal plot;
%  2. To delete an existing timespan, simply right-click on a particular bar
% Don't forget to click "Save" periodically to save your changes.
  
% MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI M-file for msprof_choose_peak_timecourses_gui.fig
%      MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI, by itself, creates a new MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI or raises the existing
%      singleton*.
%
%      H = MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI returns the handle to a new MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI or the handle to
%      the existing singleton*.
%
%      MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI.M with the given input arguments.
%
%      MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI('Property','Value',...) creates a new MSPROF_CHOOSE_PEAK_TIMECOURSES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msprof_choose_peak_timecourses_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msprof_choose_peak_timecourses_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msprof_choose_peak_timecourses_gui

% Last Modified by GUIDE v2.5 01-Jun-2010 17:38:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msprof_choose_peak_timecourses_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @msprof_choose_peak_timecourses_gui_OutputFcn, ...
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


% --- Executes just before msprof_choose_peak_timecourses_gui is made visible.
function msprof_choose_peak_timecourses_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to msprof_choose_peak_timecourses_gui (see VARARGIN)

% Choose default command line output for msprof_choose_peak_timecourses_gui
handles.output = hObject;

basename = '';
options = struct;
curarg = 1;
navigate = false;
if ~isempty(varargin)
  if isnumeric(varargin{1})
    % This is a "navigate to given m/z" call from the command line
    navigate = true;
    mz = varargin{1};
    trangeShow = [];
    curarg = 2;
    if (length(varargin) >= curarg && isnumeric(varargin{curarg}))
      trangeShow = varargin{curarg};
      curarg = curarg+1;
    end
    inputdata = getappdata(handles.figMain,'inputdata');
    if ~isempty(inputdata)
      % It's been initialized, so just navigate now
      mzMean = [inputdata.nnmfresults.mzMeanCorrected];
      [~,index] = min(abs(mzMean-mz));
      setappdata(handles.figMain,'trangeShow',trangeShow);
      mcpt_update_display(handles,index);
      return
    end
  else
    % New GUI call
    curarg = 1;
    if ischar(varargin{curarg})
      basename = varargin{curarg};
      curarg = curarg+1;
    end
  end
end
if length(varargin) >= curarg && isstruct(varargin{curarg})
  options = varargin{curarg};
end
if ~isempty(basename)
  options.basename = basename;
end
options = default(options,'basename','msprof_ecs');
options = default(options,'fit_C2',true,...
  'nnmf_results_filename',[options.basename '_nnmf.mat'],...
  'mzregistration_filename',[options.basename '_mzregistration.mat'],...
  'tregistration_filename',[options.basename '_tregistration.mat'],...
  'files_filename',[options.basename '_filelist.mat'],...
  'compounds_filename',[options.basename '_compounds.mat'],...
  'n_xaxis_pixels', 1000);

%
% Process the input data
%
inputdata = load(options.nnmf_results_filename);
n_files = length(inputdata.file);
n_peaks = length(inputdata.nnmfresults);
options = default(options,'xaxis_is_time',~isfield(inputdata,'xc'));
if options.xaxis_is_time
  options.xaxisstr = 'Time (min)';
else
  options.xaxisstr = inputdata.xaxisstr;
end
% Load info about m/z registration (needed only for "raw" display)
registration.mzShift = zeros(1,n_files);
if exist(options.mzregistration_filename,'file')
  tmp = load(options.mzregistration_filename);
  registration.mzShift = tmp.meanShiftMzI;
end
% Get information about standards, blanks, etc.
tmp = load(options.files_filename);
inputdata.standardfileindex = tmp.standardfileindex;
inputdata.blankindex = tmp.blankindex;
isotopeinfo = tmp.isotopeinfo;
inputdata.samplenames = tmp.samplenames;
% Load data about previously-defined compound spans (including
% auto-generated spans)
have_timeSpans = false;
if exist(options.compounds_filename,'file')
  cf = load(options.compounds_filename);
  if cf.randKey == inputdata.randKey
    timeSpans = cf.timeSpans;
    have_timeSpans = true;
  else
    warndlg(['If you save, file ' filenameout ' will be overwritten, because the keys do not match']);
  end
  if ~isequal(options.xaxisstr,cf.xaxisstr)
    warndlg([filenameout ' defines compounds with a different x-axis than is currently being used; if you save, file will be overwritten.']);
    have_timeSpans = false;
  end
end
if ~have_timeSpans
  timeSpans = cell(1,n_peaks);
end
% Assign a color to each sample
col = zeros(n_files,3);
for k = 1:n_files
  col(k,:) = unique_color(k+1,n_files+1); % +1 avoids using black (save for the standard)
end

% Store the data
%setappdata(hObject,'nnmf_results_filename',which(options.nnmf_results_filename));
setappdata(hObject,'options',options);
setappdata(hObject,'inputdata',inputdata);
setappdata(hObject,'timeSpans',timeSpans);
setappdata(hObject,'col',col);
setappdata(hObject,'isotopeinfo',isotopeinfo);
setappdata(hObject,'registration',registration);

% Prepare for resampling
if options.xaxis_is_time
  xrange = [max(cellfun(@(t) t(1),inputdata.tc)) min(cellfun(@(t) t(end),inputdata.tc))];
  xc = inputdata.tc;
else
  xrange = [max(cellfun(@min,inputdata.xc)),min(cellfun(@max,inputdata.xc))];
  xc = inputdata.xc;
end
xI = linspace(xrange(1),xrange(2),options.n_xaxis_pixels);
irs = zeros(n_files,length(xI));  % indices, breaks between resampled measurements
for k = 1:n_files
  keepFlag = ~isnan(xc{k});
  irs(k,:) = interp1(xc{k}(keepFlag),find(keepFlag),xI);
end
setappdata(handles.figMain,'irs',irs);
setappdata(handles.figMain,'xrange',xrange);

% Prepare isotopologue information
% This is copied from msprof_charge_isotopes_linear
n_elements = length(isotopeinfo);
Cindex = strmatch('C',{isotopeinfo.symbol},'exact');
dC = diff(isotopeinfo(Cindex).isotope_masses(1:2)); % mass increment for charge of 1
dmzc = cell(1,n_elements);
ac = cell(1,n_elements);
for i = 1:n_elements
  if (length(isotopeinfo(i).isotope_masses) > 1)
    dmzc{i} = isotopeinfo(i).isotope_masses(2:end) - ...
      isotopeinfo(i).base_mass;
    ac{i} = isotopeinfo(i).isotope_abundance(2:end) / ...
      isotopeinfo(i).isotope_abundance(1);
  end
end
l = cellfun(@length,dmzc);
isotopeIndex = make_counting_index(l);
dmz = cat(2,dmzc{:});
a = cat(2,ac{:});
% Manually append the 13C2 peak, but initialize it to zero amplitude
% (we'll fix this later)
if options.fit_C2
  dmz = [dmz 2*dC];
  a = [a 0];
  isotopeIndex = [isotopeIndex Cindex];
end
isotopologuedata = struct('dmz',dmz,'a',a,'isotopeIndex',isotopeIndex);
setappdata(hObject,'isotopologuedata',isotopologuedata);

% Prepare the controls and display
mzMean = [inputdata.nnmfresults.mzMeanCorrected];
[mzMeanS,peakRank] = sort(mzMean);
s = [num2str(mzMeanS(:)),repmat('   (',n_peaks,1),num2str(peakRank(:)),repmat(')',n_peaks,1)];
set(handles.popMz,'String',s);
setappdata(handles.figMain,'peakRank',peakRank);
axes(handles.axesSpectrum)
plot(inputdata.i2mz(1:length(inputdata.Imz_max)),inputdata.Imz_max,'k');
axis tight; set(gca,'TickDir','out')
install_mouse_event_handler(handles.axesTimeCourse,'down',@mcpt_tc_start);
install_mouse_event_handler(handles.axesTimeCourse,'up',@mcpt_tc_stop);
mcpt_update_display(handles,1);  % Start with biggest peak

% Create the legend (only need to do this once). Also put it in a place
% where it won't interfere with the parent axis
pos = get(handles.axesTimeCourse,'Position');
hleg = legend(handles.axesTimeCourse,inputdata.samplenames{:},'Location','EastOutside');
posleg = get(hleg,'Position');
posleg(1) = pos(1)+pos(3)+0.1*posleg(3);
set(hleg,'Position',posleg);

% Update handles structure
guidata(hObject, handles);

if navigate
  % This was called with "navigate" syntax, and only now are we ready to
  % execute that request
  setappdata(handles.figMain,'trangeShow',trangeShow);
  [mind,index] = min(abs(mzMean-mz));
  mcpt_update_display(handles,index);
end

% UIWAIT makes msprof_choose_peak_timecourses_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = msprof_choose_peak_timecourses_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popMz.
function popMz_Callback(hObject, eventdata, handles)
% hObject    handle to popMz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popMz contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popMz
  mzIndex = get(hObject,'Value');
  peakRank = getappdata(handles.figMain,'peakRank');
  index = peakRank(mzIndex);
  mcpt_update_display(handles,index);
  

% --- Executes during object creation, after setting all properties.
function popMz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editRank_Callback(hObject, eventdata, handles)
% hObject    handle to editRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRank as text
%        str2double(get(hObject,'String')) returns contents of editRank as a double
  peakRank = round(str2double(get(hObject,'String')));
  mcpt_update_display(handles,peakRank);
  


% --- Executes during object creation, after setting all properties.
function editRank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnPrevious.
function btnPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  peakRank = round(str2double(get(handles.editRank,'String')));
  mcpt_update_display(handles,peakRank-1);


% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% hObject    handle to btnNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  peakRank = round(str2double(get(handles.editRank,'String')));
  mcpt_update_display(handles,peakRank+1);


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  options = getappdata(handles.figMain,'options');
  timeSpans = getappdata(handles.figMain,'timeSpans');
  inputdata = getappdata(handles.figMain,'inputdata');
  randKey = inputdata.randKey;
  xaxisstr = options.xaxisstr;
  save(options.compounds_filename,'timeSpans','randKey','xaxisstr');

function mcpt_update_display(handles,index)
  inputdata = getappdata(handles.figMain,'inputdata');
  col = getappdata(handles.figMain,'col');
  peakRank = getappdata(handles.figMain,'peakRank');
  thispeak = inputdata.nnmfresults(index);
  xrange = getappdata(handles.figMain,'xrange');
  options = getappdata(handles.figMain,'options');
  irs = getappdata(handles.figMain,'irs');
  state = {'off','on'};
  set(handles.btnPrevious,'Enable',state{(index > 1)+1});
  set(handles.btnNext,'Enable',state{(index < length(inputdata.nnmfresults))+1});
  
  % Update the spectrum axis
  hline = findobj(handles.axesSpectrum,'Tag','thisMz');
  if ~isempty(hline)
    delete(hline);
  end
  rng = thispeak.Imzrange;
  rng = rng(1):rng(2);
  y = max(inputdata.Imz_max(rng));
  line([0 0]+thispeak.mzMeanCorrected,[0 y],'Parent',handles.axesSpectrum,'Color','r','Tag','thisMz');
  % Update the m/z pulldown
  pdindex = find(peakRank == index);
  set(handles.popMz,'Value',pdindex);
  % Update the rank edit box
  set(handles.editRank,'String',num2str(index));
  % Update the fractional intensity box
  fI = thispeak.totIntensity / inputdata.nnmfresults(1).totIntensity;
  set(handles.textTotI,'String',sprintf('Fractional intensity: %g',fI));
  % Update the charge box
  set(handles.editCharge,'String',num2str(thispeak.charge))
  % Plot the time course
  delete(get(handles.axesTimeCourse,'Children'))
  n_files = length(thispeak.iProf);
  x = linspace(xrange(1),xrange(2),size(irs,2)-1);
  for k = 1:n_files
    % Interpolate in a way that preserves the sum of counts
    Irs = diff(interp1(cumsum(thispeak.iProf{k}),irs(k,:)));
    if get(handles.checkboxLogScale,'Value')
      Irs = sqrt(Irs);
    end
    line(x,Irs,'Parent',handles.axesTimeCourse,'Color',col(k,:),'HitTest','off');
  end
  set(handles.axesTimeCourse,'XTick',[],'XLim',xrange);
  set(handles.axesTimeCourseLabel,'String',options.xaxisstr);
  mcpt_update_timespans(handles,index)
  mcpt_update_raw(handles,index)
  
function mcpt_update_timespans(handles,index)
  timeSpans = getappdata(handles.figMain,'timeSpans');
  thisTS = timeSpans{index};
  xrange = getappdata(handles.figMain,'xrange');
  set(handles.axesTimeSpans,'XLim',xrange);
  % Delete any current lines in this axis
  delete(get(handles.axesTimeSpans,'Children'))
  % Plot any spans
  if ~isempty(thisTS)
    hline = line(thisTS,[1;1]*(1:size(thisTS,2)),'Parent',handles.axesTimeSpans,'Color','k','ButtonDownFcn',@mcpt_delete_tS);
    set(handles.axesTimeSpans,'YLim',[0 size(thisTS,2)]+0.5);
    % Determine whether we should highlight one in red (this would mean
    % we are in mz & span callback mode)
    trangeShow = getappdata(handles.figMain,'trangeShow');
    if ~isempty(trangeShow)
      redFlag = all(thisTS == repmat(trangeShow(:),1,size(thisTS,2)),1);
      set(hline(redFlag),'Color','r');
    end
    setappdata(handles.figMain,'trangeShow',[]); % clear the selection
  end

function mcpt_update_raw(handles,index)
  hax = [handles.axesRaw handles.axesRawp1 handles.axesRawp2];
  for i = 1:length(hax)
    cla(hax(i));
  end
  if get(handles.checkboxShowRaw,'Value') == 0
    return
  end
  logscale = get(handles.checkboxLogScale,'Value');
  raw = getappdata(handles.figMain,'raw');
  fileIndex = get(handles.popupmenuFile,'Value');
  inputdata = getappdata(handles.figMain,'inputdata');
  isotopologuedata = getappdata(handles.figMain,'isotopologuedata');
  xrange = getappdata(handles.figMain,'xrange');
  irs = getappdata(handles.figMain,'irs');
  mzBase = inputdata.nnmfresults(index).mzMeanCorrected;
  charge = round(str2double(get(handles.editCharge,'String')));
  if (charge < 1)
    charge = 1;
  end
  set(handles.editCharge,'String',num2str(charge));
  mzList = mzBase + isotopologuedata.dmz/charge;
  mzIndex = inputdata.mz2i([mzBase mzList]) - inputdata.standardMzIShift;
  width = 20;
  [regions,regionI] = split_into_contiguous_regions(round(mzIndex),width);
  n_regions = size(regions,2);
  biggest_dmz = zeros(1,n_regions);
  a = [1 isotopologuedata.a];
  dmz = [0 isotopologuedata.dmz];
  for regionIndex = 1:n_regions
    thisindx = regionI == regionIndex;
    thisa = a(thisindx);
    thisdmz = dmz(thisindx);
    [~,maxaIndex] = max(thisa);
    biggest_dmz(regionIndex) = thisdmz(maxaIndex);
  end
  for axIndex = 1:length(hax)
    rng = regions(1,axIndex):regions(2,axIndex);
    S = raw{fileIndex}(rng);
    yl = inputdata.i2mz(rng+inputdata.standardMzIShift);
    [i,j,s] = ssparse_find(S);
    img = accumarray([i(:) j(:)],s(:));
    l = length(inputdata.tc{fileIndex});
    img(1,end+1:l) = 0;  % fill out to full size, if necessary
    img1 = diff(interp1(cumsum(img,2)',irs(fileIndex,:)))';
    if logscale
      img1 = sqrt(img1);
    end
    if (axIndex == 1)
      cl = max(img1(:));
    end
    image('Parent',hax(axIndex),'CData',img1,'CDataMapping','scaled','YData',yl,'XData',xrange);
    set(hax(axIndex),'XLim',xrange,'YLim',yl([1 end]),'CLim',[0 cl],'XTick',[]);
    line('Parent',hax(axIndex),'XData',xrange,'YData',(mzBase+biggest_dmz(axIndex))*[1 1],'Color','k','LineStyle','--');
  end
  


function handled = mcpt_tc_start(sender,event)
  cp = get(sender,'CurrentPoint');
  line([0 0]+cp(1,1),[0 0]+cp(1,2),'Parent',sender,'Color','k','LineWidth',2,'EraseMode','xor','Tag','dragline');
  handles = guidata(sender);
  install_mouse_event_handler(handles.axesTimeCourse,'move',@mcpt_tc_drag);
  handled = true;

function handled = mcpt_tc_drag(sender,event)
  cp = get(sender,'CurrentPoint');
  hline = findobj(sender,'Tag','dragline');
  if isempty(hline)
    return
  end
  xd = get(hline,'XData');
  set(hline,'XData',[xd(1) cp(1,1)]);
  handled = true;
  
function handled = mcpt_tc_stop(sender,event)
  hline = findobj(sender,'Tag','dragline');
  xd = get(hline,'XData');
  if iscell(xd)
    % Something weird happened...
    return
  end
  handles = guidata(sender);
  index = round(str2double(get(handles.editRank,'String')));
  timeSpans = getappdata(handles.figMain,'timeSpans');
  thisTS = timeSpans{index};
  if isempty(thisTS)
    thisTS = xd(:);
  else
    thisTS = [thisTS xd(:)];
  end
  timeSpans{index} = thisTS;
  setappdata(handles.figMain,'timeSpans',timeSpans);
  delete(hline)
  uninstall_mouse_event_handler(handles.axesTimeCourse,'move',@mcpt_tc_drag);
  mcpt_update_timespans(handles,index)
  handled = true;

function mcpt_delete_tS(sender,eventdata)
  handles = guidata(sender);
  timeSpans = getappdata(handles.figMain,'timeSpans');
  index = round(str2double(get(handles.editRank,'String')));
  thisTS = timeSpans{index};
  xd = get(sender,'XData');
  tSindex = thisTS(1,:) == xd(1) & thisTS(2,:) == xd(2);
  thisTS(:,tSindex) = [];
  if isempty(thisTS)
    thisTS = zeros(2,0);
  end
  timeSpans{index} = thisTS;
  setappdata(handles.figMain,'timeSpans',timeSpans);
  mcpt_update_timespans(handles,index)


% --- Executes on button press in checkboxAuto.
function checkboxAuto_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAuto


% --- Executes on button press in btnFormula.
function btnFormula_Callback(hObject, eventdata, handles)
% hObject    handle to btnFormula (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  index = round(str2double(get(handles.editRank,'String')));
  inputdata = getappdata(handles.figMain,'inputdata');
  isotopeinfo = getappdata(handles.figMain,'isotopeinfo');
  thispeak = inputdata.nnmfresults(index);
  formula_gui(thispeak.mzMeanCorrected*thispeak.charge,isotopeinfo,thispeak.nI,thispeak.nIcov);
  


% --- Executes on button press in checkboxLogScale.
function checkboxLogScale_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLogScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLogScale
  peakRank = round(str2double(get(handles.editRank,'String')));
  mcpt_update_display(handles,peakRank);


% --- Executes on button press in checkboxShowRaw.
function checkboxShowRaw_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowRaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowRaw
raw = getappdata(handles.figMain,'raw');
if isempty(raw)
  % This is the first call, must load the file data
  inputdata = getappdata(handles.figMain,'inputdata');
  n_files = length(inputdata.file);
  s = getappdata(handles.figMain,'registration');
  s = default(s,'mzShift',zeros(1,n_files),'tShift',zeros(1,n_files));
  inputdata = getappdata(handles.figMain,'inputdata');
  mzShift = -s.mzShift;
  mzShift = mzShift - min(mzShift);
  [Mc,mz2i,i2mz,tc] = msprof_load_set(inputdata.file,struct('mzi_shift',mzShift,'t_shift',-s.tShift));
  Mp = cell(size(Mc));
  fprintf('Now transposing inputs...');
  for i = 1:length(Mc)
    Mp{i} = ssparse_swap(Mc{i});
  end
  disp('..done');
  setappdata(handles.figMain,'raw',Mp);
  set(handles.popupmenuFile,'String',inputdata.samplenames,'Value',1);
end
peakRank = round(str2double(get(handles.editRank,'String')));
mcpt_update_display(handles,peakRank);



% --- Executes on selection change in popupmenuFile.
function popupmenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuFile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFile
  peakRank = round(str2double(get(handles.editRank,'String')));
  mcpt_update_display(handles,peakRank);


% --- Executes during object creation, after setting all properties.
function popupmenuFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCharge_Callback(hObject, eventdata, handles)
% hObject    handle to editCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCharge as text
%        str2double(get(hObject,'String')) returns contents of editCharge as a double
   peakRank = round(str2double(get(handles.editRank,'String')));
   mcpt_update_raw(handles,peakRank)

% --- Executes during object creation, after setting all properties.
function editCharge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
