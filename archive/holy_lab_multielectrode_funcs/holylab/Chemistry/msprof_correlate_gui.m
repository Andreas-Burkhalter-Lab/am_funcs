function varargout = msprof_correlate_gui(varargin)
% MSPROF_CORRELATE_GUI: correlate neuronal responses to ligand abundance
%
% This is a GUI for exploring the relationship between neuronal responses
% and ligand abundance.
%
% Syntax:
%   msprof_correlate_gui(drdata,concdata,lsqerr)
% where
%   drdata and concdata are as described in msprof_monotonic_dr;
%   lsqerr is the output of msprof_correlate_dr.
%
% Usage:
% When the GUI starts, the top panel shows responses of cells to the
% different stimuli.  Each cell is in a column.
%   If you click on a cell, then you see in the middle panel a plot of how
% well each compound can account for the responses of the cell. The
% horizontal axis represents the different compounds, and the vertical axis
% represents the error in fitting the cell's response to that compound's
% concentration profile. The most interesting compounds are the ones with
% smallest error.
%   If you then click on a particular compound, you can see a plot of the
% cell's firing rate to all of the stimuli (at all concentrations), sorted
% in order of increasing concentration of the ligand across all tested
% samples.  The firing rate should look quite monotonic, or have errors
% (e.g., spike-amplitude decrement) that are well-understood. Also note
% that in this window, click-drag zooms in on a range of compounds;
% right-clicking restores the full view.
%   As an alternative to clicking on a cell in the top panel to get things
% started, you can select a particular compound's m/z in the upper right
% pulldown menu. Then the middle panel plots the error as a function of
% cell #. However, note that since different cells have different degrees
% of noise, it is not straightforward to interpret the resulting plot.
%   Finally, the small "MS" button in the lower right allows you to inspect
% the original MS data on the chosen compound---this can be useful in
% determining whether you're "in the noise" or not.
%
% See also: MSPROF_CORRELATE_DR, MSPROF_MONOTONIC_DR.

% Copyright 2009-2010 by Timothy E. Holy

% MSPROF_CORRELATE_GUI M-file for msprof_correlate_gui.fig
%      MSPROF_CORRELATE_GUI, by itself, creates a new MSPROF_CORRELATE_GUI or raises the existing
%      singleton*.
%
%      H = MSPROF_CORRELATE_GUI returns the handle to a new MSPROF_CORRELATE_GUI or the handle to
%      the existing singleton*.
%
%      MSPROF_CORRELATE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSPROF_CORRELATE_GUI.M with the given input arguments.
%
%      MSPROF_CORRELATE_GUI('Property','Value',...) creates a new MSPROF_CORRELATE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msprof_correlate_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msprof_correlate_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msprof_correlate_gui

% Last Modified by GUIDE v2.5 03-Feb-2011 04:26:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msprof_correlate_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @msprof_correlate_gui_OutputFcn, ...
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


% --- Executes just before msprof_correlate_gui is made visible.
function msprof_correlate_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to msprof_correlate_gui (see VARARGIN)

% Choose default command line output for msprof_correlate_gui
handles.output = hObject;

drdata = varargin{1};
concdata = varargin{2};
lsqerr = varargin{3};

% Sort the compounds in order of most abundant to least abundant
maxconc = max(concdata.I,[],1);
[maxconcs,sortOrder] = sort(maxconc,'descend');
concdata.mz = concdata.mz(sortOrder);
concdata.xrange = concdata.xrange(:,sortOrder);
concdata.I = concdata.I(:,sortOrder);
if isfield(concdata,'nI')
  concdata.nI = concdata.nI(:,sortOrder);
end
lsqerr = lsqerr(:,sortOrder);
% Prepare the popup, by also sorting in order of m/z
[mzs,peakRank] = sort(concdata.mz);
n_peaks = length(mzs);
s = [num2str(mzs(:)),repmat('   [',n_peaks,1),num2str(concdata.xrange(1,peakRank)'),repmat(',',n_peaks,1),num2str(concdata.xrange(2,peakRank)'),repmat(']   (',n_peaks,1),num2str(peakRank(:)),repmat(')',n_peaks,1)];
set(handles.popMz,'String',s);
setappdata(handles.figMain,'peakRank',peakRank);

% Assign colors to each stimulus
if isfield(drdata,'uniquecolors')
  col = drdata.uniquecolors;
else
  n_stimuli = length(drdata.uniquelabels);
  col = distinguishable_colors(n_stimuli);
%   col = zeros(n_stimuli,3);
%   for i = 1:n_stimuli
%     col(i,:) = unique_color(i,n_stimuli);
%   end
end
if ~isfield(drdata,'uniquelinestyles')
  drdata.uniquelinestyles = repmat({'-'},1,length(drdata.uniquelabels));
end

% Store the data in the figure
setappdata(handles.figMain,'drdata',drdata);
setappdata(handles.figMain,'concdata',concdata);
setappdata(handles.figMain,'lsqerr',lsqerr);
setappdata(handles.figMain,'col',col);

% Get & save isotopeinfo, if applicable
if isfield(concdata,'nI')
  isotopeinfo = isotopes;
  Sindx = strmatch('S',{isotopeinfo.symbol},'exact');
  setappdata(handles.figMain,'Sindx',Sindx);
  % Create a second axis for plotting # sulfurs
  hax = copyobj(handles.axesAllCorrelations,handles.figMain);
  set(hax,'HitTest','off','Color','none','YAxisLocation','right');
  handles.axesAllCorrelations2 = hax;
else
  set(handles.checkboxSulfur,'Visible','off');
end

% Plot the all-cells data
% mx = max(abs(drdata.relconc(:)));
%mdexplore_im(drdata.relconc,struct('plotfunc',@mspcg_pickcell,'hax',handles.axesAllCells,'clim',mx*[-1 1],'ylabel',{drdata.uniquelabels}));
mdexplore_im(drdata.relconc,struct('plotfunc',@mspcg_pickcell,'hax',handles.axesAllCells,'ylabel',{drdata.uniquelabels}));
% mdexplore_im(bsxfun(@rdivide,drdata.relconc,max(drdata.relconc,[],1)),struct('plotfunc',@mspcg_pickcell,'hax',handles.axesAllCells,'ylabel',{drdata.uniquelabels})); % normalized
cm = colormap_pm(drdata.relconc);
colormap(handles.axesAllCells,cm);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes msprof_correlate_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = msprof_correlate_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function mspcg_pickcell(cellIndex,row)
  handles = guidata(gcbf);
  setappdata(handles.figMain,'thisCell',cellIndex);
  % Plot the correlation to individual compounds
  lsqerr = getappdata(handles.figMain,'lsqerr');
  hl = plot(handles.axesAllCorrelations,lsqerr(cellIndex,:),'b');
  set(hl,'HitTest','off','Tag','lsqerr');
  set(handles.axesAllCorrelations,'Box','off','TickDir','out')
  if isappdata(handles.figMain,'Sindx')
    Sindx = getappdata(handles.figMain,'Sindx');
    concdata = getappdata(handles.figMain,'concdata');
    hs = plot(handles.axesAllCorrelations2,concdata.nI(Sindx,:),'r');
    visible = get(handles.checkboxSulfur,'Value');
    offon = {'off','on'};
    set(hs,'HitTest','off','Visible',offon{1+visible})
    set(handles.axesAllCorrelations2,'Color','none','YLim',[-0.4 2],'YAxisLocation','right','Box','off','TickDir','out','HitTest','off','Visible',offon{1+visible})
  end
  set(handles.axesAllCorrelations,'ButtonDownFcn',@(src,event) mspcg_pickcorrelation(src,event,'compoundIndex'));
  cla(handles.axesThisCellSorted);
  cla(handles.axesThisCellConc);
  set(handles.textSingleLigand,'String','')
  set(handles.textCorrelation,'String',sprintf('Cell %d correlation',cellIndex));
  setappdata(handles.figMain,'compoundIndex',[]);
  % Plot the raw responses of this cell
  mspcg_plotCell(handles)
    
  
function mspcg_pickcorrelation(src,event,tag)
  % This function handles mouse button clicks within the correlations
  % window. Left-click: choose compound; left-click-drag: zoom in;
  % right-click (or any other click): zoom out to full size
  % This picks a compound based on the correlations to a cell, or picks a
  % cell based on correlations to a compound
  handles = guidata(gcbf);
  selType = get(gcbf,'SelectionType');
  if (selType(1) == 'n')
    % Left-click: determine whether this is a zoom or a selection
    startp = get(handles.axesAllCorrelations,'CurrentPoint');
    startp = startp(1,1:2);
    rect = rbbox;
    endp = get(handles.axesAllCorrelations,'CurrentPoint');
    endp = endp(1,1:2);
    if (rect(3) < 3)
      % This is a click. Choose a compound.
      [pos,hl,selIndex] = findpoint(startp,handles.axesAllCorrelations);
      setappdata(handles.figMain,tag,selIndex);
      if strcmp(tag,'thisCell')
        mspcg_plotCell(handles)
      end
      mspcg_updateCellLigand(handles)
    else
      % This is a drag. Zoom in.
      hax = handles.axesAllCorrelations;
      if isfield(handles,'axesAllCorrelations2')
        hax(2) = handles.axesAllCorrelations2;
      end
      set(hax,'XLim',sort([startp(1) endp(1)]));
    end
  else
    % Right-click: zoom out to full view
    hl = findobj(handles.axesAllCorrelations,'Tag','lsqerr');
    if ishandle(hl)
      xd = get(hl,'XData');
      hax = handles.axesAllCorrelations;
      if isfield(handles,'axesAllCorrelations2')
        hax(2) = handles.axesAllCorrelations2;
      end
      set(hax,'XLim',[1 length(xd)]);
    end
  end
  
function mspcg_updateCellLigand(handles)
  cellIndex = getappdata(handles.figMain,'thisCell');
  compoundIndex = getappdata(handles.figMain,'compoundIndex');
  drdata = getappdata(handles.figMain,'drdata');
  concdata = getappdata(handles.figMain,'concdata');
  %  This is code from msprof_correlate_dr
  % Find the common stimuli
  clabel = intersect(concdata.stimtag,drdata.uniquelabels);
  flagdr = ismember(drdata.stimlabel,clabel);
  flagc = ismember(concdata.stimtag,clabel);
  % For each (matching) entry in drdata, find the corresponding entry in concdata
  mapc2dr = findainb(drdata.stimlabel(flagdr),concdata.stimtag(flagc));
  % Compute the concentration of this ligand in each presented stimulus
  concfac = drdata.stimconc(flagdr);
  tmpI0 = concdata.I(mapc2dr,compoundIndex);
  tmpI = tmpI0 .* concfac(:); % Contains estimate of ligand concentration in each presented sample
  [sI,sortOrder] = sort(tmpI);  % Sort them into increasing concentration
  % Plot the responses in order, sorted by concentration
  tmpr = drdata.rates(flagdr,cellIndex);
  tmpre = drdata.rateerrs(flagdr,cellIndex);
  errorbar(handles.axesThisCellSorted,1:length(tmpr),tmpr(sortOrder),tmpre(sortOrder));
  set(handles.axesThisCellSorted,'XLim',[0 length(tmpr)+1]);
  % Also plot against the actual concentration
  cla(handles.axesThisCellConc);
  hold(handles.axesThisCellConc,'on');
  flagudr = ismember(drdata.uniquelabels,clabel);
  for i = 1:length(flagudr)
    if ~flagudr(i)
      continue
    end
    rindx = strmatch(drdata.uniquelabels{i},drdata.stimlabel,'exact');
    cindx = strmatch(drdata.uniquelabels{i},concdata.stimtag,'exact');
    c = concdata.I(cindx,compoundIndex) * drdata.stimconc(rindx);
    hl = errorbar(handles.axesThisCellConc,c,drdata.rates(rindx,cellIndex),...
      drdata.rateerrs(rindx,cellIndex));
    set(hl,'Color',drdata.uniquecolors(i,:),'LineStyle',drdata.uniquelinestyles{i});
  end
%   errorbar(handles.axesThisCellConc,sI,tmpr(sortOrder),tmpre(sortOrder));
  xl = get(handles.axesThisCellConc,'XLim');
  xl(2) = sI(end);
%   xl(1) = max(xl(1),1e-5*sI(end));
  xl(1) = 1e-5*sI(end);
  xl(2) = xl(2)*1.1;
  set(handles.axesThisCellConc,'XLim',xl,'XScale','log')
  set(handles.textSingleLigand,'String',sprintf('Cell %d, m/z %.4f, xrange [%.3g  %.3g]',cellIndex,concdata.mz(compoundIndex),concdata.xrange(:,compoundIndex)));
  
function mspcg_plotCell(handles)
  drdata = getappdata(handles.figMain,'drdata');
  cellIndex = getappdata(handles.figMain,'thisCell');
  col = getappdata(handles.figMain,'col');
  % Plot the responses as a function of concentration
  [ul,tmp,stimI] = unique(drdata.stimlabel);
  cl = agglabel(stimI);
  cla(handles.axesThisCell);
  hold(handles.axesThisCell,'on')
  hl = zeros(1,length(cl));
  for sI = 1:length(cl)
    hl(sI) = errorbar(drdata.stimconc(cl{sI}),drdata.rates(cl{sI},cellIndex),drdata.rateerrs(cl{sI},cellIndex),'Parent',handles.axesThisCell);
    set(hl(sI),'Color',col(sI,:),'LineStyle',drdata.uniquelinestyles{sI});%,'LineWidth',lw(sI));
  end
  set(handles.textThisCell,'String',sprintf('Cell %d',cellIndex));
  legend(hl,drdata.uniquelabels,'Location','EastOutside');
  xl = [min(drdata.stimconc) max(drdata.stimconc)];
  xl = log10(xl);
  xl = xl + [-1 1]*0.05*diff(xl);
  set(handles.axesThisCell,'XScale','log','XLim',10.^xl);

  
  


% --- Executes on button press in pushbuttonMS.
function pushbuttonMS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  concdata = getappdata(handles.figMain,'concdata');
  compoundIndex = getappdata(handles.figMain,'compoundIndex');
  msprof_choose_peak_timecourses_gui(concdata.mz(compoundIndex),concdata.xrange(:,compoundIndex));
  


% --- Executes on selection change in popMz.
function popMz_Callback(hObject, eventdata, handles)
% hObject    handle to popMz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popMz contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popMz
  selectionIndex = get(hObject,'Value');
  if isfield(handles,'axesAllCorrelations2')
    set(handles.checkboxSulfur,'Value',0);
    set(handles.axesAllCorrelations2,'Visible','off')
  end
  peakRank = getappdata(handles.figMain,'peakRank');
  compoundIndex = peakRank(selectionIndex);
  setappdata(handles.figMain,'compoundIndex',compoundIndex);
  % Draw the lsqerr for all cells, given this compound
  lsqerr = getappdata(handles.figMain,'lsqerr');
  drdata = getappdata(handles.figMain,'drdata');
  mse = mean(drdata.relconcerr.^2,1);
  lsqerrnorm = lsqerr(:,compoundIndex) ./ sqrt(mse(:));  % normalize by each cells RMSE
  %hl = plot(handles.axesAllCorrelations,lsqerr(:,compoundIndex));
  hl = plot(handles.axesAllCorrelations,lsqerrnorm);
  set(hl,'HitTest','off');
  set(handles.axesAllCorrelations,'YScale','log','ButtonDownFcn',@(src,event) mspcg_pickcorrelation(src,event,'thisCell'));
  cla(handles.axesThisCellSorted);
  concdata = getappdata(handles.figMain,'concdata');
  set(handles.textCorrelation,'String',sprintf('Compound %g [%g %g] correlation',concdata.mz(compoundIndex),concdata.xrange(:,compoundIndex)));
  setappdata(handles.figMain,'thisCell',[]);

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


% --- Executes on button press in checkboxSulfur.
function checkboxSulfur_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSulfur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSulfur
  offon = {'off','on'};
  value = get(hObject,'Value');
  set(handles.axesAllCorrelations2,'Visible',offon{1+value});
  hc = get(handles.axesAllCorrelations2,'Children');
  set(hc,'Visible',offon{1+value});
