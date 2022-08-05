function varargout = mschoosetrials_gui(varargin)
% MSCHOOSETRIALS_GUI: GUI for selecting valid mass spectra
% Syntax:
%   trials = mschoosetrials_gui(msdata)
%   trials = mschoosetrials_gui(msdata,m)
% where
%   msdata is a mass spec structure of the type loaded by MSLOAD;
%   m (optional) is the full raw mass spectrum for this fraction;
% and
%   trials is a vector of valid trial numbers (can be empty).
%
% See also: MSLOAD.
  
% Copyright 2005 by Timothy E. Holy
  
% MSCHOOSETRIALS_GUI M-file for mschoosetrials_gui.fig
%      MSCHOOSETRIALS_GUI, by itself, creates a new MSCHOOSETRIALS_GUI or raises the existing
%      singleton*.
%
%      H = MSCHOOSETRIALS_GUI returns the handle to a new MSCHOOSETRIALS_GUI or the handle to
%      the existing singleton*.
%
%      MSCHOOSETRIALS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSCHOOSETRIALS_GUI.M with the given input arguments.
%
%      MSCHOOSETRIALS_GUI('Property','Value',...) creates a new MSCHOOSETRIALS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mschoosetrials_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mschoosetrials_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help mschoosetrials_gui

% Last Modified by GUIDE v2.5 11-Feb-2005 09:51:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mschoosetrials_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mschoosetrials_gui_OutputFcn, ...
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


% --- Executes just before mschoosetrials_gui is made visible.
function mschoosetrials_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mschoosetrials_gui (see VARARGIN)

% Choose default command line output for mschoosetrials_gui
%handles.output = hObject;

figure(hObject);
% The following is preparation for a callback in which ms spectra are
% shown on a per-trial basis if the user clicks a trial marker
% Compute spectra of identified peaks for each trial
msdata = varargin{1};
npeaks = length(msdata.peaks);
spec = [];
for i = 1:npeaks
  ntrials = length(msdata.peaks(i).trial);
  for j = 1:ntrials
    ctrial = msdata.peaks(i).trial(j);
    if (ctrial > length(spec))
      spec(ctrial).mz = [];
      spec(ctrial).mag = [];
    end
    spec(ctrial).mz(end+1) = msdata.peaks(i).mz(j);
    spec(ctrial).mag(end+1) = msdata.peaks(i).mag(j);
  end
end
% If the user supplied raw msdata, parse that as well
if (length(varargin) > 1)
  m = varargin{2};
  dm = diff(m(:,2));
  ibreak = [1;find(dm < 0)+1;length(dm)+2];
  ntrials = length(ibreak)-1;
  for i = 1:ntrials
    spec(i).rawmz = m(ibreak(i):ibreak(i+1)-1,:);
  end
end
options.ncorr = 10;

% Find the largest peaks (since we want more than one of these for
% correlation purposes, mschoosepeak is not the best choice)
mxmag = zeros(1,npeaks);
for i = 1:npeaks
  mxmag(i) = max(msdata.peaks(i).mag);
end
[smxmag,mxmagindx] = sort(mxmag,'descend');
mxmagindx = mxmagindx(1:min(options.ncorr,npeaks));
bj = mxmagindx(1);

set(handles.donebutton,'UserData',msdata.peaks(bj).trial);

%
% Draw the GUI
%
% Plot the amplitude vs. trial number
axes(handles.axes1);
[xtmp,ytmp] = addbreaks(msdata.peaks(bj).trial,msdata.peaks(bj).mag);
h = plot(xtmp,ytmp,'b');
set(h,'HitTest','off');
hold on
% Also plot as discrete points, so can set up callback for each trial
ntrials = length(msdata.peaks(bj).trial);
for i = 1:ntrials
  ctrial = msdata.peaks(bj).trial(i);  
  hpoint = plot(ctrial,msdata.peaks(bj).mag(i),'bo');
  spec(ctrial).trialnum = ctrial;
  set(hpoint,'UserData',spec(ctrial),'MarkerFaceColor','b',...
    'MarkerSize',12,...
    'ButtonDownFcn',@msct_showspectrum);
end
  
% Draw the two range selection lines
ymax = 1.1*max(msdata.peaks(bj).mag);
xmax = max(msdata.peaks(bj).trial);
hline = line([0.5 xmax+0.5; 0.5 xmax+0.5],[0 0; ymax ymax],'Color','k', ...
             'LineStyle','-','EraseMode','xor');
drag_line(hline(1),struct('type','v'));
drag_line(hline(2),struct('type','v'));
xlimbuf = ceil(0.05*xmax);

set(gca,'XLim',[0-xlimbuf xmax+1+xlimbuf],'YLim',[0 ymax]);


% Compute the correlation in relative magnitude of the largest peaks
peakmag = nan(length(mxmagindx),ntrials);  % Each column is one trial; each
                                           % peak is one row
for i = 1:size(peakmag,1)
  peakmag(i,msdata.peaks(mxmagindx(i)).trial) = msdata.peaks(mxmagindx(i)).mag';
end
[rtrial,ptrial] = corrcoef(peakmag);
axes(handles.corrax);
pos = get(handles.corrax,'Position');
imagesc(rtrial);%,[0 1]);
hcb = colorbar('SouthOutside');


% Title the plot
set(handles.filenametext,'String',msdata.filename);

% Update handles structure
handles.rangelines = hline;
guidata(hObject, handles);

% UIWAIT makes mschoosetrials_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mschoosetrials_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
if isempty(handles)
  varargout{1} = -1;
  return;
end
trials = [];
if (get(handles.discardbox,'Value') == get(handles.discardbox,'Min'))
  xdata = get(handles.rangelines,'XData');
  xdata = cat(2,xdata{:});
  trials = get(handles.donebutton,'UserData');
  trials(find(trials < min(xdata))) = [];
  trials(find(trials > max(xdata))) = [];
end
varargout{1} = trials;
close(handles.figure1);


% --- Executes on button press in discardbox.
function discardbox_Callback(hObject, eventdata, handles)
% hObject    handle to discardbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of discardbox


% --- Executes on button press in donebutton.
function donebutton_Callback(hObject, eventdata, handles)
% hObject    handle to donebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


function msct_showspectrum(hObject,event)
  spec = get(hObject,'UserData');
  figure;
  hold on
  if isfield(spec,'rawmz')
    stem(spec.rawmz(:,2),spec.rawmz(:,1),'k.');
  end
  stem(spec.mz,spec.mag,'r.');
  xlabel('m/z');
  ylabel('Ion current');
  set(gca,'TickDir','out')
  title(['Trial ' num2str(spec.trialnum)]);
